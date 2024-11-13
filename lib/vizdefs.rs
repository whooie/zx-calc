use tabbycat::attributes::Color;

pub const FONT: &str = "DejaVu Sans";
pub const FONTSIZE: f64 = 10.0; // pt
pub const SQUARE_HEIGHT: f64 = 0.15; // in
pub const CIRCLE_HEIGHT: f64 = 0.20; // in
pub const NODE_MARGIN: f64 = 0.025; // in

// pub const Z_COLOR: Color = Color::White;
// pub const X_COLOR: Color = Color::Gray50;
// pub const H_COLOR: Color = Color::White;

pub const Z_COLOR  : Color = Color::Rgb(115, 150, 250); // blue
pub const X_COLOR  : Color = Color::Rgb(230, 115, 125); // red
pub const H_COLOR  : Color = Color::Rgb(250, 205, 115); // yellow
pub const H_WIRE   : Color = Color::Rgb(68,  136, 255); // lighter blue
pub const STAR_WIRE: Color = Color::Rgb(255, 136, 68 ); // orange

