//! The structure gathering id of sketched genome and its probminhash


/// The structure to be embedded in Hnsw.
/// The distance declared in Hnsw will be be on IdSketch, operating on the field id
/// TODO : possibly replace the String id by a hash to spare memory in Hnsw
pub struct IdSketch{
    /// id of genome Sketched
    id : String
    /// probminhash signature
    sig : Vec<u32>

} // end of IdSketch