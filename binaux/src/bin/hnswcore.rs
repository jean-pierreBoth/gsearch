//! This binary is dedicated to coreset computations on data stored in Hnsw created by crate [hnsw_rs](https://crates.io/crates/hnsw_rs).  
//! It is cloned from crate [fromhnsw](https://crates.io/crates/fromhnsw) in workspace [coreset](https://crates.io/crates/coreset)
//!
//! 1. The simple default command just computes a coreset
//! The command is:
//! **hnscore  --dir (-d) dirname  --fname (-f) hnswname  --typename (-t) typename**.
//!
//!  The following arguments are mandatory:
//! - dirname : directory where hnsw files reside
//! - hnswname : name used for naming the 2 hnsw related files: name.hnsw.data and name.hnsw.graph
//! - typename : can be u16, u32, u64, f32, f64, i16, i32, i64 depending on the Distance type.
//!
//! At the end of coreset computations all data are re-scanned and dispatched to  nearest facility center.  
//!  A csv file named *corest.csv* is dumped in current directory.  
//!     Each line consists in 2 DataId, the first one identifies a data point and the second the DataId of the center of its corresponding facility.
//!
//! 2. The clustercore subcommand parses the following parameters optional parameters:
//!
//! - cluster :
//!     number of cluster asked in the Kmedoid pass. This argument is optional and defaults to 0. If non zero, it asks for a Kmedoid end pass on the facilities created by
//!     the coreset algorithm, generating the *nbcluster* asked for.
//!     At the end of coreset/cluster computations all data are re-scanned and dispatched to nearest cluster center.  
//!   
//!     A csv file named *clustercorest.csv* is dumped in current directory.  
//!     Each line consists in 2 DataId, the first one identifies a data point and the second the DataId of the center of its corresponding cluster.  
//!
//! - fraction :
//!      The size of the coreset generated will be around fraction * size of data.
//!      A point can be sampled many times, in this case the sampled points are merged and their weight added.
//!      So the fraction should be set to a value slightly superior to the one desired.  
//!      For large data, it can/must be set to a lower value.  
//!
//! The following optional arguments are related to the first Bmor pass of the algorithm. They are explained in detail in Bmor documentation):  
//!  - beta:  defaults to 2. The accepted cost evolves as beta^iter during iterations. Increasing beta makes the accepted cost greater, and sp
//!     reduces the number of facilities generated.  
//!  - gamma: defaults to 2. Increasing gamma allocates a greater number of facilites.
//!
//! command is :  
//! **hnscore  --dir (-d) dirname  --fname (-f) hnswname  --typename (-t) typename  clustercore --cluster nbcluster [--beta (-b)  value] [--gamma (-g) value]**
//!
//!
//!
//! Note: It is easy to add any adhoc type T  by adding a line in [get_datamap()].  
//! The type T used in hnsw must satisfy: 'static + Clone + Sized + Send + Sync + std::fmt::Debug

use cpu_time::ProcessTime;
use std::time::{Duration, SystemTime};

use coreset::prelude::*;

use clap::{Arg, ArgAction, ArgMatches, Command};
use std::default::Default;

use hnsw_rs::datamap::*;
use hnsw_rs::prelude::*;

use fromhnsw::getdatamap::get_typed_datamap;
use fromhnsw::hnswiter::HnswMakeIter;

use quantiles::ckms::CKMS;
use rand::Rng;

//========================================
// Parameters

#[derive(Debug, Clone)]
struct HnswParams {
    dir: String,
    hname: String,
    // type name of base Vector in Hnsw to reload
    tname: String,
}

impl HnswParams {
    pub fn new(hdir: &str, hname: &str, tname: &str) -> Self {
        HnswParams {
            dir: hdir.to_string(),
            hname: hname.to_string(),
            tname: tname.to_string(),
        }
    }
}

//
/// Coreset parameters
#[derive(Copy, Clone, Debug)]
pub struct CoresetParams {
    beta: f32,
    gamma: f32,
    freduc: f32,
    // if clusterization is askerd after coreset, nbcluster contains the number of cluster required, oterwise is set to 0
    nbcluster: usize,
}

impl CoresetParams {
    #[allow(unused)]
    fn new(beta: f32, gamma: f32, freduc: f32, nbcluster: usize) -> CoresetParams {
        CoresetParams {
            beta,
            gamma,
            freduc,
            nbcluster,
        }
    }
    //
    fn get_beta(&self) -> f32 {
        self.beta
    }

    //
    fn get_gamma(&self) -> f32 {
        self.gamma
    }
    //
    fn get_reduction(&self) -> f32 {
        self.freduc
    }
    // return number of cluster asked for if any , else 0
    fn get_cluster(&self) -> usize {
        self.nbcluster
    }
} // end of impl CoresetParams

impl Default for CoresetParams {
    fn default() -> Self {
        CoresetParams {
            beta: 2.,
            gamma: 2.,
            freduc: 0.11,
            nbcluster: 0,
        }
    }
}

//===========================================================

fn parse_coreset_cmd(matches: &ArgMatches) -> Result<CoresetParams, anyhow::Error> {
    log::debug!("in  parse_coreset_cmd");
    //
    let gamma = *matches.get_one::<f32>("gamma").unwrap();
    let beta = *matches.get_one::<f32>("beta").unwrap();
    log::info!("gamma : {}", gamma);
    let params = CoresetParams::new(
        beta,
        gamma,
        *matches.get_one::<f32>("reduction").unwrap(),
        *matches.get_one::<usize>("cluster").unwrap(),
    );
    //
    log::info!("got CoresetParams : {:?}", params);
    //
    Ok(params)
}

//============================================================================================

/// This function dispatch its call to get_typed_datamap::\<T\> according to type T
/// The cuurent function dispatch to u16, u32, u64, i32, i64, f32 and f64 according to typename.
/// For another type, the functio is easily modifiable.  
/// The only constraints on T comes from hnsw and is T: 'static + Clone + Sized + Send + Sync + std::fmt::Debug
pub fn get_datamap(directory: String, basename: String, typename: &str) -> anyhow::Result<DataMap> {
    //
    let datamap = match typename {
        "u16" => get_typed_datamap::<u16>(directory, basename),
        "u32" => get_typed_datamap::<u32>(directory, basename),
        "u64" => get_typed_datamap::<u64>(directory, basename),
        "f32" => get_typed_datamap::<f32>(directory, basename),
        "f64" => get_typed_datamap::<f64>(directory, basename),
        "i32" => get_typed_datamap::<i32>(directory, basename),
        "i64" => get_typed_datamap::<i64>(directory, basename),
        _ => {
            log::error!(
                "get_datamap : unimplemented type, type received : {}",
                typename
            );
            std::panic!("get_datamap : unimplemented type");
        }
    };
    //
    log::info!("returning DataMap for type : {}", typename);
    //
    datamap
}

//===========================================================

#[allow(unused)]
fn check_distances<T, Dist>(distance: &Dist, datamap: &DataMap)
where
    T: Send + Sync + Clone + std::fmt::Debug,
    Dist: Distance<T> + Sync + Send + Clone + Default,
{
    log::info!("check_distances");
    let mut q_dist = CKMS::<f32>::new(0.01);
    let nb_to_load = 1000.min(datamap.get_nb_data());
    let mut data_v = Vec::<Vec<T>>::with_capacity(nb_to_load);
    let mut keys = datamap.get_dataid_iter();
    while data_v.len() < nb_to_load {
        let dataid = keys.next().unwrap();
        let v = datamap.get_data::<T>(dataid).unwrap();
        data_v.push(Vec::from(v));
    }
    // now sample indexes and compute dists
    let nb_sample = 50000;
    let mut sampled = 0;
    let mut rng = rand::rng();
    while sampled < nb_sample {
        let i = rng.random_range(0..nb_to_load);
        let j = rng.random_range(0..nb_to_load);
        if i != j {
            let dist = distance.eval(&data_v[i], &data_v[j]);
            q_dist.insert(dist);
            sampled += 1;
        }
    }
    //
    println!("statistics on sampled distances : ");
    println!(
        "\n distance quantiles at  0.005 : {:.2e} , 0.01 :  {:.2e} , 0.025 : {:.2e}, 0.25 : {:.2e}, 0.5 : {:.2e}, 0.75 : {:.2e}   0.99 : {:.2e}\n",
        q_dist.query(0.005).unwrap().1,
        q_dist.query(0.01).unwrap().1,
        q_dist.query(0.025).unwrap().1,
        q_dist.query(0.25).unwrap().1,
        q_dist.query(0.5).unwrap().1,
        q_dist.query(0.75).unwrap().1,
        q_dist.query(0.99).unwrap().1
    );
}

//===========================================================

fn coreset1<T, Dist>(coreparams: &CoresetParams, datamap: &DataMap) -> usize
where
    T: Send + Sync + Clone + std::fmt::Debug,
    Dist: Distance<T> + Sync + Send + Clone + Default,
{
    //
    let dist_name = std::any::type_name::<Dist>();
    log::info!("coreset1 instantiated for distance = {:?}", dist_name);
    let distance = Dist::default();
    //
    check_distances(&distance, datamap);
    //
    println!("\n\n entering coreset + our kmedoids");
    println!("==================================");
    let cpu_start = ProcessTime::now();
    let sys_now = SystemTime::now();
    //
    let beta = coreparams.get_beta().into();
    let gamma = coreparams.get_gamma().into();
    let freduc: f64 = coreparams.get_reduction().into();
    let nb_data = datamap.get_nb_data();
    //
    let iter_producer = HnswMakeIter::<T>::new(datamap);
    // now do we have only coreset computation or also clusterization
    if coreparams.get_cluster() == 0 {
        let k = 10; // as we have 10 classes, but this gives a lower bound
        let mut core1 = Coreset1::<usize, T, Dist>::new(k, nb_data, beta, gamma, distance.clone());
        //
        let res = core1.make_coreset(&iter_producer, freduc);
        if res.is_err() {
            log::error!("construction of coreset1 failed");
        }
        let coreset = res.unwrap();
        let dump_res = coreset.dump();
        if dump_res.is_err() {
            log::error!("error occurred dump coreset in csv file");
        }
    } else {
        // we must do coreset + clusterization
        let bmor_arg = BmorArg::new(nb_data, beta, gamma);
        let nb_max_kmedoid_iter = 15;
        let mut clustercoreset =
            ClusterCoreset::<usize, T>::new(coreparams.get_cluster(), freduc, bmor_arg);
        clustercoreset.compute(distance.clone(), nb_max_kmedoid_iter, &iter_producer);
        clustercoreset.dispatch(&distance, &iter_producer);
    }
    // dump a csv with membership.
    //
    let cpu_time: Duration = cpu_start.elapsed();
    println!(
        "  sys time(ms) {:?} cpu time(ms) {:?}",
        sys_now.elapsed().unwrap().as_millis(),
        cpu_time.as_millis()
    );
    1
} // end of

//===========================================================

fn main() {
    //
    println!("\n ************** initializing logger *****************\n");
    env_logger::Builder::from_default_env().init();
    log::info!("logger initialized from default environment");
    //
    let core_params: CoresetParams;

    //
    let coresetcmd = Command::new("coreset")
        .about("Coreset approximate clustering")
        .version("0.3.2")
        .arg(
            Arg::new("beta")
                .required(false)
                .short('b')
                .long("beta")
                .default_value("2.0")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(f32))
                .help("beta"),
        )
        .arg(
            Arg::new("gamma")
                .required(false)
                .short('g')
                .long("gamma")
                .default_value("2.0")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(f32))
                .help("gamma"),
        )
        .arg(
            Arg::new("reduction")
                .required(false)
                .short('f')
                .long("fraction")
                .default_value("0.1")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(f32))
                .help("reduction"),
        )
        .arg(
            Arg::new("cluster")
                .required(false)
                .short('c')
                .long("cluster")
                .default_value("0")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(usize))
                .help("cluster"),
        );
    //
    // global command
    // =============
    //
    let matches = Command::new("hnswcore")
        .about("Approximate Clustering via Streaming Coreset Construction and HNSW")
        .arg_required_else_help(true)
        .arg(
            Arg::new("dir")
                .long("dir")
                .short('d')
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(String))
                .required(true)
                .help("expecting a directory name"),
        )
        .arg(
            Arg::new("fname")
                .long("fname")
                .short('f')
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(String))
                .required(true)
                .help("expecting a file  basename"),
        )
        .arg(
            Arg::new("typename")
                .short('t')
                .long("type")
                .value_parser(clap::value_parser!(String))
                .required(true)
                .help("expecting a directory name"),
        )
        .subcommand(coresetcmd)
        .get_matches();
    //
    // retrieve HnswPathParams
    //
    let hdir = matches
        .get_one::<String>("dir")
        .expect("dir argument needed");
    let hname = matches
        .get_one::<String>("fname")
        .expect("hnsw base name needed");
    let tname: &String = matches
        .get_one::<String>("typename")
        .expect("typename required");
    //
    let hparams = HnswParams::new(hdir, hname, tname);
    log::info!("received parameters: {:?}", hparams);
    //
    // parse coreset parameters
    //
    if let Some(core_match) = matches.subcommand_matches("coreset") {
        log::debug!("subcommand for coreset parameters");
        let res = parse_coreset_cmd(core_match);
        match res {
            Ok(params) => {
                core_params = params;
            }
            _ => {
                log::error!("parsing coreset command failed");
                println!("exiting with error {}", res.err().as_ref().unwrap());
                //  log::error!("exiting with error {}", res.err().unwrap());
                std::process::exit(1);
            }
        }
    } else {
        core_params = CoresetParams::default();
    }
    log::debug!("coreset params : {:?}", core_params);
    //
    // Datamap Creation
    //
    let typename = &hparams.tname;
    let datamap = get_datamap(hparams.dir, hparams.hname, typename);
    if datamap.is_err() {
        log::error!(
            "datamap could not be constructed : {}",
            datamap.err().unwrap()
        );
        std::process::exit(1);
    }
    let datamap = datamap.unwrap();
    // Distance instanciation
    let distname = &datamap.get_distname();
    let typename = &datamap.get_data_typename();
    let short_name = distname.split("::").last().unwrap();
    //
    let _res = if short_name == "DistL1" {
        match typename.as_str() {
            "u16" => coreset1::<u16, DistL1>(&core_params, &datamap),
            "u32" => coreset1::<u32, DistL1>(&core_params, &datamap),
            "f32" => coreset1::<f32, DistL1>(&core_params, &datamap),
            _ => panic!("not yet"),
        }
    } else if short_name == "DistL2" {
        match typename.as_str() {
            "u16" => coreset1::<u16, DistL2>(&core_params, &datamap),
            "u32" => coreset1::<u32, DistL2>(&core_params, &datamap),
            "f32" => coreset1::<f32, DistL2>(&core_params, &datamap),
            _ => panic!("not yet implemented type"),
        }
    } else if short_name == "DistCosine" {
        match typename.as_str() {
            "f32" => coreset1::<f32, DistCosine>(&core_params, &datamap),
            "f64" => coreset1::<f64, DistCosine>(&core_params, &datamap),
            _ => panic!("not yet implemented type"),
        }
    } else if short_name == "DistHamming" {
        match typename.as_str() {
            "u16" => coreset1::<u16, DistHamming>(&core_params, &datamap),
            "u32" => coreset1::<u32, DistHamming>(&core_params, &datamap),
            "f32" => coreset1::<f32, DistHamming>(&core_params, &datamap),
            _ => panic!("not yet implemented type"),
        }
    } else {
        log::error!("not yet implemented distance,  distance : {:?}", distname);
        panic!("not yet implemented distance");
    };

    // need a macro or switch function to dispatch on types
}
