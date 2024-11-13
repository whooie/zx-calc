use crate::phase::Phase;

/// A single spider.
#[derive(Copy, Clone, Debug)]
pub enum Spider {
    /// A Z-spider, parameterized by a real phase.
    Z(Phase),
    /// An X-spider, parameterized by a real phase.
    X(Phase),
}

impl Spider {
    /// Create a new, phaseless Z-spider.
    pub fn z() -> Self { Self::Z(Phase::zero()) }

    /// Create a new Z-spider with π phase.
    pub fn z_pi() -> Self { Self::Z(Phase::pi()) }

    /// Create a new Z-spider with π/2 phase.
    pub fn z_pi2() -> Self { Self::Z(Phase::pi2()) }

    /// Create a new Z-spider with π/4 phase.
    pub fn z_pi4() -> Self { Self::Z(Phase::pi4()) }

    /// Create a new Z-spider with π/8 phase.
    pub fn z_pi8() -> Self { Self::Z(Phase::pi8()) }

    /// Create a new Z-spider with phase `(a / b) × 2π`.
    pub fn z_frac(a: i64, b: i64) -> Self { Self::Z(Phase::new(a, b)) }

    /// Create a new, phaseless X-spider.
    pub fn x() -> Self { Self::X(Phase::zero()) }

    /// Create a new X-spider with π phase.
    pub fn x_pi() -> Self { Self::X(Phase::pi()) }

    /// Create a new X-spider with π/2 phase.
    pub fn x_pi2() -> Self { Self::X(Phase::pi2()) }

    /// Create a new X-spider with π/4 phase.
    pub fn x_pi4() -> Self { Self::X(Phase::pi4()) }

    /// Create a new X-spider with π/8 phase.
    pub fn x_pi8() -> Self { Self::X(Phase::pi8()) }

    /// Cre ate a new X-spider with phase `(a / b) × 2π`.
    pub fn x_frac(a: i64, b: i64) -> Self { Self::X(Phase::new(a, b)) }

    /// Return `true` if `self` is `Z`.
    pub fn is_z(&self) -> bool { matches!(self, Self::Z(_)) }

    /// Return `true` if `self` is `Z` and the phase satisfies some predicate.
    pub fn is_z_and<F>(&self, pred: F) -> bool
    where F: FnOnce(Phase) -> bool
    {
        match self {
            Self::Z(ph) => pred(*ph),
            _ => false,
        }
    }

    /// Return `true` if `self` is `X`.
    pub fn is_x(&self) -> bool { matches!(self, Self::X(_)) }

    /// Return `true` if `self` is `X` and the phase satisfies some predicate.
    pub fn is_x_and<F>(&self, pred: F) -> bool
    where F: FnOnce(Phase) -> bool
    {
        match self {
            Self::X(ph) => pred(*ph),
            _ => false,
        }
    }

    /// Return `true` if the phase of the spider is 0.
    pub fn has_defarg(&self) -> bool {
        match self {
            Self::Z(ph) | Self::X(ph) => *ph == Phase::zero(),
        }
    }

    /// Return the associated phase.
    pub fn phase(&self) -> Phase {
        match self {
            Self::Z(ph) | Self::X(ph) => *ph,
        }
    }

    /// Return `true` if `self` is a spider with the given phase, modulo 2π.
    pub fn has_phase(&self, phase: Phase) -> bool {
        match self {
            Self::Z(ph) | Self::X(ph) => *ph == phase,
        }
    }

    /// Return `true` if `self` and `other` are both spiders with the same
    /// color.
    pub fn is_same_color(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (Self::Z(_), Self::Z(_)) | (Self::X(_), Self::X(_))
        )
    }

    /// Return `true` if `self` and `other` are both spiders with the same color
    /// and their phases satisfy a predicate.
    pub fn is_same_color_and<F>(&self, other: &Self, pred: F) -> bool
    where F: FnOnce(Phase, Phase) -> bool
    {
        match (self, other) {
            (Self::Z(ph1), Self::Z(ph2)) | (Self::X(ph1), Self::X(ph2))
                => pred(*ph1, *ph2),
            _ => false,
        }
    }

    /// Return `true` if `self` and `other` are both spiders with different
    /// color.
    pub fn is_diff_color(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (Self::Z(_), Self::X(_)) | (Self::X(_), Self::Z(_))
        )
    }

    /// Return `true` if `self` and `other` are both spiders with different
    /// color and their phases satisfy a predicate.
    pub fn is_diff_color_and<F>(&self, other: &Self, pred: F) -> bool
    where F: FnOnce(Phase, Phase) -> bool
    {
        match (self, other) {
            (Self::Z(ph1), Self::X(ph2)) | (Self::X(ph1), Self::Z(ph2))
                => pred(*ph1, *ph2),
            _ => false,
        }
    }

}

