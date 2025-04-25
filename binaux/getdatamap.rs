//! This module gets a DataMap from hnsw dump

use log;

use std::path::Path;

use anyhow;

use hnsw_rs::datamap::*;

/// reloads a datamap and checks for type T.
pub fn get_typed_datamap<T: 'static + std::fmt::Debug>(
    directory: String,
    basename: String,
) -> anyhow::Result<DataMap> {
    //
    let path = Path::new(&directory);
    let res = DataMap::from_hnswdump::<T>(path, &basename);
    if res.is_err() {
        log::error!(
            "get_datamap, could not get datamap from hnsw, directory {}, basename : {}",
            directory,
            basename
        );
    }
    let datamap = res.unwrap();
    let t_name = datamap.get_data_typename();
    // check type
    let check_type = datamap.check_data_type::<T>();
    if !check_type {
        log::error!(
            "bad type name. registered type name : {}, you asked for {}",
            t_name,
            std::any::type_name::<T>().to_string()
        )
    }
    //
    Ok(datamap)
}
