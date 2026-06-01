//! Minimal BinaryCIF decoder: MessagePack container + the 7 column encodings.
//! Spec: https://github.com/molstar/BinaryCIF/blob/master/encoding.md
use anyhow::{anyhow, bail, Result};
use rmpv::Value;

enum Num {
    I(Vec<i64>),
    F(Vec<f64>),
}

fn as_i(n: Num) -> Result<Vec<i64>> {
    match n {
        Num::I(v) => Ok(v),
        Num::F(_) => bail!("expected integer array in encoding chain"),
    }
}

fn field<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    v.as_map()?
        .iter()
        .find(|(k, _)| k.as_str() == Some(key))
        .map(|(_, val)| val)
}

fn kind_of(enc: &Value) -> Option<&str> {
    field(enc, "kind").and_then(Value::as_str)
}

/// Interpret a raw little-endian byte buffer as a typed number array (ByteArray).
fn byte_array(data: &[u8], ty: i64) -> Result<Num> {
    let n = |sz: usize| data.len() / sz;
    Ok(match ty {
        1 => Num::I(data.iter().map(|&b| b as i8 as i64).collect()),
        2 => Num::I((0..n(2)).map(|i| i16::from_le_bytes([data[2 * i], data[2 * i + 1]]) as i64).collect()),
        3 => Num::I((0..n(4)).map(|i| i32::from_le_bytes(data[4 * i..4 * i + 4].try_into().unwrap()) as i64).collect()),
        4 => Num::I(data.iter().map(|&b| b as i64).collect()),
        5 => Num::I((0..n(2)).map(|i| u16::from_le_bytes([data[2 * i], data[2 * i + 1]]) as i64).collect()),
        6 => Num::I((0..n(4)).map(|i| u32::from_le_bytes(data[4 * i..4 * i + 4].try_into().unwrap()) as i64).collect()),
        32 => Num::F((0..n(4)).map(|i| f32::from_le_bytes(data[4 * i..4 * i + 4].try_into().unwrap()) as f64).collect()),
        33 => Num::F((0..n(8)).map(|i| f64::from_le_bytes(data[8 * i..8 * i + 8].try_into().unwrap())).collect()),
        other => bail!("unsupported ByteArray type {other}"),
    })
}

/// Accumulate packed values: values equal to the type's limit are summed with the
/// following value(s) until a non-limit value ends the run.
fn integer_packing_accumulate(raw: &[i64], byte_count: i64, unsigned: bool) -> Vec<i64> {
    let (upper, lower): (i64, i64) = if unsigned {
        (if byte_count == 1 { 0xFF } else { 0xFFFF }, 0)
    } else if byte_count == 1 {
        (0x7F, -0x80)
    } else {
        (0x7FFF, -0x8000)
    };
    let mut out = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        let mut value = 0i64;
        let mut t = raw[i];
        while (unsigned && t == upper) || (!unsigned && (t == upper || t == lower)) {
            value += t;
            i += 1;
            t = raw[i];
        }
        out.push(value + t);
        i += 1;
    }
    out
}

/// Decode IntegerPacking directly from a raw byte buffer (when it is innermost).
fn integer_packing(data: &[u8], byte_count: i64, unsigned: bool) -> Vec<i64> {
    let raw: Vec<i64> = match (byte_count, unsigned) {
        (1, true) => data.iter().map(|&b| b as i64).collect(),
        (1, false) => data.iter().map(|&b| b as i8 as i64).collect(),
        (2, true) => (0..data.len() / 2).map(|i| u16::from_le_bytes([data[2 * i], data[2 * i + 1]]) as i64).collect(),
        _ => (0..data.len() / 2).map(|i| i16::from_le_bytes([data[2 * i], data[2 * i + 1]]) as i64).collect(),
    };
    integer_packing_accumulate(&raw, byte_count, unsigned)
}

/// The innermost encoding consumes the raw byte buffer.
fn consume_bytes(data: &[u8], enc: &Value) -> Result<Num> {
    match kind_of(enc) {
        Some("ByteArray") => byte_array(data, field(enc, "type").and_then(Value::as_i64).ok_or_else(|| anyhow!("ByteArray.type"))?),
        Some("IntegerPacking") => {
            let byte_count = field(enc, "byteCount").and_then(Value::as_i64).ok_or_else(|| anyhow!("byteCount"))?;
            let unsigned = field(enc, "isUnsigned").and_then(Value::as_bool).unwrap_or(false);
            Ok(Num::I(integer_packing(data, byte_count, unsigned)))
        }
        other => bail!("unexpected innermost encoding {other:?}"),
    }
}

/// Subsequent encodings transform an already-decoded number array.
fn apply(enc: &Value, num: Num) -> Result<Num> {
    Ok(match kind_of(enc) {
        Some("IntegerPacking") => {
            let byte_count = field(enc, "byteCount").and_then(Value::as_i64).ok_or_else(|| anyhow!("byteCount"))?;
            let unsigned = field(enc, "isUnsigned").and_then(Value::as_bool).unwrap_or(false);
            Num::I(integer_packing_accumulate(&as_i(num)?, byte_count, unsigned))
        }
        Some("RunLength") => {
            let v = as_i(num)?;
            let mut out = Vec::new();
            for pair in v.chunks_exact(2) {
                for _ in 0..pair[1] {
                    out.push(pair[0]);
                }
            }
            Num::I(out)
        }
        Some("Delta") => {
            let origin = field(enc, "origin").and_then(Value::as_i64).unwrap_or(0);
            let v = as_i(num)?;
            let mut val = origin;
            Num::I(v.into_iter().map(|d| {
                val += d;
                val
            }).collect())
        }
        Some("FixedPoint") => {
            let factor = field(enc, "factor").and_then(Value::as_f64).ok_or_else(|| anyhow!("factor"))?;
            Num::F(as_i(num)?.into_iter().map(|d| d as f64 / factor).collect())
        }
        Some("IntervalQuantization") => {
            let min = field(enc, "min").and_then(Value::as_f64).ok_or_else(|| anyhow!("min"))?;
            let max = field(enc, "max").and_then(Value::as_f64).ok_or_else(|| anyhow!("max"))?;
            let steps = field(enc, "numSteps").and_then(Value::as_i64).ok_or_else(|| anyhow!("numSteps"))?;
            let delta = (max - min) / (steps - 1) as f64;
            Num::F(as_i(num)?.into_iter().map(|d| min + delta * d as f64).collect())
        }
        other => bail!("unexpected transform encoding {other:?}"),
    })
}

fn decode_numeric(data: &[u8], encodings: &[Value]) -> Result<Num> {
    let mut rev = encodings.iter().rev();
    let first = rev.next().ok_or_else(|| anyhow!("empty encoding chain"))?;
    let mut num = consume_bytes(data, first)?;
    for enc in rev {
        num = apply(enc, num)?;
    }
    Ok(num)
}

fn string_array(col_data: &[u8], enc: &Value) -> Result<Vec<Option<String>>> {
    let data_enc = field(enc, "dataEncoding").and_then(Value::as_array).ok_or_else(|| anyhow!("dataEncoding"))?;
    let offset_enc = field(enc, "offsetEncoding").and_then(Value::as_array).ok_or_else(|| anyhow!("offsetEncoding"))?;
    let offset_bytes = field(enc, "offsets").and_then(Value::as_slice).ok_or_else(|| anyhow!("offsets"))?;
    let string_data = field(enc, "stringData").and_then(Value::as_str).ok_or_else(|| anyhow!("stringData"))?;

    let offsets = as_i(decode_numeric(offset_bytes, offset_enc)?)?;
    let indices = as_i(decode_numeric(col_data, data_enc)?)?;
    let bytes = string_data.as_bytes();
    let substrings: Vec<&str> = (0..offsets.len().saturating_sub(1))
        .map(|i| std::str::from_utf8(&bytes[offsets[i] as usize..offsets[i + 1] as usize]).unwrap_or(""))
        .collect();

    Ok(indices
        .into_iter()
        .map(|ix| if ix < 0 { None } else { substrings.get(ix as usize).map(|s| s.to_string()) })
        .collect())
}

/// A decoded column. `None` entries are masked (CIF `.`/`?`).
pub enum ColumnData {
    Str(Vec<Option<String>>),
    Int(Vec<Option<i64>>),
    Float(Vec<Option<f64>>),
}

impl ColumnData {
    fn apply_mask(&mut self, mask: &[i64]) {
        match self {
            ColumnData::Str(v) => v.iter_mut().zip(mask).for_each(|(x, &m)| if m != 0 { *x = None }),
            ColumnData::Int(v) => v.iter_mut().zip(mask).for_each(|(x, &m)| if m != 0 { *x = None }),
            ColumnData::Float(v) => v.iter_mut().zip(mask).for_each(|(x, &m)| if m != 0 { *x = None }),
        }
    }

    pub fn as_str(&self, row: usize) -> Option<String> {
        match self {
            ColumnData::Str(v) => v.get(row).cloned().flatten(),
            ColumnData::Int(v) => v.get(row).copied().flatten().map(|x| x.to_string()),
            ColumnData::Float(v) => v.get(row).copied().flatten().map(|x| x.to_string()),
        }
    }

    pub fn as_f64(&self, row: usize) -> Option<f64> {
        match self {
            ColumnData::Float(v) => v.get(row).copied().flatten(),
            ColumnData::Int(v) => v.get(row).copied().flatten().map(|x| x as f64),
            ColumnData::Str(v) => v.get(row).cloned().flatten().and_then(|s| s.trim().parse().ok()),
        }
    }

    pub fn as_i64(&self, row: usize) -> Option<i64> {
        match self {
            ColumnData::Int(v) => v.get(row).copied().flatten(),
            ColumnData::Float(v) => v.get(row).copied().flatten().map(|x| x as i64),
            ColumnData::Str(v) => v.get(row).cloned().flatten().and_then(|s| s.trim().parse().ok()),
        }
    }
}

pub struct Category {
    pub name: String,
    pub row_count: usize,
    columns: Vec<(String, ColumnData)>,
}

impl Category {
    pub fn column(&self, name: &str) -> Option<&ColumnData> {
        self.columns.iter().find(|(n, _)| n == name).map(|(_, c)| c)
    }
}

fn decode_column(col: &Value) -> Result<(String, ColumnData)> {
    let name = field(col, "name").and_then(Value::as_str).ok_or_else(|| anyhow!("column.name"))?.to_string();
    let data = field(col, "data").ok_or_else(|| anyhow!("column.data"))?;
    let data_bytes = field(data, "data").and_then(Value::as_slice).ok_or_else(|| anyhow!("data.data"))?;
    let encs = field(data, "encoding").and_then(Value::as_array).ok_or_else(|| anyhow!("data.encoding"))?;

    let mut decoded = if encs.first().map(kind_of) == Some(Some("StringArray")) {
        ColumnData::Str(string_array(data_bytes, &encs[0])?)
    } else {
        match decode_numeric(data_bytes, encs)? {
            Num::I(v) => ColumnData::Int(v.into_iter().map(Some).collect()),
            Num::F(v) => ColumnData::Float(v.into_iter().map(Some).collect()),
        }
    };

    if let Some(mask_v) = field(col, "mask") {
        if !mask_v.is_nil() {
            let md = field(mask_v, "data").and_then(Value::as_slice).ok_or_else(|| anyhow!("mask.data"))?;
            let me = field(mask_v, "encoding").and_then(Value::as_array).ok_or_else(|| anyhow!("mask.encoding"))?;
            decoded.apply_mask(&as_i(decode_numeric(md, me)?)?);
        }
    }
    Ok((name, decoded))
}

/// Decode the first data block of a BinaryCIF buffer into its categories.
pub fn decode_first_block(bytes: &[u8]) -> Result<Vec<Category>> {
    let root = rmpv::decode::read_value(&mut &bytes[..]).map_err(|e| anyhow!("messagepack: {e}"))?;
    let blocks = field(&root, "dataBlocks").and_then(Value::as_array).ok_or_else(|| anyhow!("dataBlocks"))?;
    let block = blocks.first().ok_or_else(|| anyhow!("no data blocks"))?;
    let categories = field(block, "categories").and_then(Value::as_array).ok_or_else(|| anyhow!("categories"))?;

    let mut out = Vec::new();
    for cat in categories {
        let name = field(cat, "name").and_then(Value::as_str).unwrap_or("").trim_start_matches('_').to_string();
        let row_count = field(cat, "rowCount").and_then(Value::as_i64).unwrap_or(0) as usize;
        let cols = field(cat, "columns").and_then(Value::as_array).ok_or_else(|| anyhow!("columns"))?;
        let columns = cols.iter().map(decode_column).collect::<Result<Vec<_>>>()?;
        out.push(Category { name, row_count, columns });
    }
    Ok(out)
}
