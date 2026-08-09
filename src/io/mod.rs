// Decoder selection is always available: `ThreadConfig` carries the user's
// preference and the CLIs consult it whether or not the parallel decoder is
// compiled in.
pub mod broker;
pub mod calibrate;
pub mod fastx;
pub mod map_info;
pub mod policy;
pub mod rad;
pub mod threads;
