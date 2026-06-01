use anyhow::Result;
use arrow::array::{ArrayRef, Float64Array, Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use serde_json::Value;
use std::collections::BTreeSet;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use structscope_features::FeatureRecord;

#[derive(Clone, Copy, PartialEq)]
enum ColKind {
    Int,
    Float,
    Str,
}

/// Infer a column type from the union of values seen for `key`.
/// All-integer -> Int64, all-numeric -> Float64, anything else -> Utf8 (JSON).
fn classify(records: &[FeatureRecord], key: &str) -> ColKind {
    let mut kind = ColKind::Int;
    for record in records {
        match record.features.get(key) {
            None | Some(Value::Null) => {}
            Some(Value::Number(n)) => {
                if !(n.is_i64() || n.is_u64()) && kind == ColKind::Int {
                    kind = ColKind::Float;
                }
            }
            Some(_) => return ColKind::Str,
        }
    }
    kind
}

fn value_to_string(value: &Value) -> String {
    match value {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

/// Write feature records as a columnar Parquet file with a stable, sorted schema.
pub fn write_parquet(path: &Path, records: &[FeatureRecord]) -> Result<()> {
    let mut keys: BTreeSet<&str> = BTreeSet::new();
    for record in records {
        for key in record.features.keys() {
            keys.insert(key.as_str());
        }
    }

    let mut fields = vec![
        Field::new("structure_id", DataType::Utf8, false),
        Field::new("source_path", DataType::Utf8, true),
    ];
    let mut columns: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(
            records.iter().map(|r| Some(r.structure_id.clone())).collect::<Vec<_>>(),
        )),
        Arc::new(StringArray::from(
            records.iter().map(|r| r.source_path.clone()).collect::<Vec<_>>(),
        )),
    ];

    for key in &keys {
        let (data_type, array): (DataType, ArrayRef) = match classify(records, key) {
            ColKind::Int => (
                DataType::Int64,
                Arc::new(Int64Array::from(
                    records.iter().map(|r| r.features.get(*key).and_then(Value::as_i64)).collect::<Vec<_>>(),
                )),
            ),
            ColKind::Float => (
                DataType::Float64,
                Arc::new(Float64Array::from(
                    records.iter().map(|r| r.features.get(*key).and_then(Value::as_f64)).collect::<Vec<_>>(),
                )),
            ),
            ColKind::Str => (
                DataType::Utf8,
                Arc::new(StringArray::from(
                    records.iter().map(|r| r.features.get(*key).map(value_to_string)).collect::<Vec<_>>(),
                )),
            ),
        };
        fields.push(Field::new(*key, data_type, true));
        columns.push(array);
    }

    let schema = Arc::new(Schema::new(fields));
    let batch = RecordBatch::try_new(schema.clone(), columns)?;
    let mut writer = ArrowWriter::try_new(File::create(path)?, schema, None)?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{Array, Float64Array, Int64Array, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use serde_json::json;
    use structscope_features::FeatureRecord;

    fn record(id: &str, features: serde_json::Value) -> FeatureRecord {
        FeatureRecord {
            structure_id: id.to_string(),
            source_path: Some(format!("/data/{id}.pdb")),
            features: features.as_object().unwrap().clone(),
        }
    }

    #[test]
    fn writes_parquet_with_inferred_types() {
        let records = vec![
            record("a", json!({ "atom_count": 10, "rg": 1.5, "centroid": [0.0, 1.0, 2.0] })),
            record("b", json!({ "atom_count": 20, "rg": 2.5, "centroid": [3.0, 4.0, 5.0] })),
        ];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("f.parquet");
        write_parquet(&path, &records).unwrap();

        let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(&path).unwrap())
            .unwrap()
            .build()
            .unwrap();
        let batch = reader.into_iter().next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);

        let col = |name: &str| batch.column(batch.schema().index_of(name).unwrap()).clone();
        let atoms = col("atom_count");
        assert_eq!(atoms.as_any().downcast_ref::<Int64Array>().unwrap().value(0), 10);
        let rg = col("rg");
        assert_eq!(rg.as_any().downcast_ref::<Float64Array>().unwrap().value(1), 2.5);
        // arrays serialize to JSON strings
        let centroid = col("centroid");
        assert_eq!(
            centroid.as_any().downcast_ref::<StringArray>().unwrap().value(0),
            "[0.0,1.0,2.0]"
        );
        let ids = col("structure_id");
        assert_eq!(ids.as_any().downcast_ref::<StringArray>().unwrap().value(0), "a");
    }
}
