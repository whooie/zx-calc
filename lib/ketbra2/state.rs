use std::{ collections::VecDeque, fmt };
use num_complex::Complex64 as C64;
use itertools::Itertools;

/// Like [`Ord`], but without the expectation that the ordering is based on the
/// whole value; this trait should produce orderings based *only* on quantum
/// state information.
pub trait StateOrd {
    fn state_cmp(&self, other: &Self) -> std::cmp::Ordering;
}

/// The state on a single wire.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum State {
    /// ∣0⟩ = ∣+⟩ + ∣–⟩
    Zero,
    /// ∣1⟩ = ∣+⟩ – ∣–⟩
    One,
    /// ∣+⟩ = ∣0⟩ + ∣1⟩
    Plus,
    /// ∣–⟩ = ∣0⟩ - ∣1⟩
    Minus,
}

impl State {
    /// Calculate the dot product `⟨self|rhs⟩`.
    pub fn dot(&self, other: &Self) -> C64 {
        use std::f64::consts::FRAC_1_SQRT_2;
        match (self, other) {
            (Self::Zero,  Self::Zero ) => 1.0.into(),
            (Self::Zero,  Self::One  ) => 0.0.into(),
            (Self::Zero,  Self::Plus ) => FRAC_1_SQRT_2.into(),
            (Self::Zero,  Self::Minus) => FRAC_1_SQRT_2.into(),
            (Self::One,   Self::Zero ) => 0.0.into(),
            (Self::One,   Self::One  ) => 1.0.into(),
            (Self::One,   Self::Plus ) => FRAC_1_SQRT_2.into(),
            (Self::One,   Self::Minus) => (-FRAC_1_SQRT_2).into(),
            (Self::Plus,  Self::Zero ) => FRAC_1_SQRT_2.into(),
            (Self::Plus,  Self::One  ) => FRAC_1_SQRT_2.into(),
            (Self::Plus,  Self::Plus ) => 1.0.into(),
            (Self::Plus,  Self::Minus) => 0.0.into(),
            (Self::Minus, Self::Zero ) => FRAC_1_SQRT_2.into(),
            (Self::Minus, Self::One  ) => (-FRAC_1_SQRT_2).into(),
            (Self::Minus, Self::Plus ) => 0.0.into(),
            (Self::Minus, Self::Minus) => 1.0.into(),
        }
    }

    /// Return `true` if `self` is a member of `basis`.
    pub fn is_basis(&self, basis: Basis) -> bool {
        matches!(
            (self, basis),
            (Self::Zero,  Basis::Z)
            | (Self::One,   Basis::Z)
            | (Self::Plus,  Basis::X)
            | (Self::Minus, Basis::X)
        )
    }

    /// Return `self` decomposed in `basis` if `self` is not already a member of
    /// `basis`.
    pub fn decomp(&self, basis: Basis) -> Option<((C64, Self), (C64, Self))> {
        if self.is_basis(basis) {
            None
        } else {
            let up = basis.up();
            let down = basis.down();
           Some( ( (self.dot(&up), up), (self.dot(&down), down) ) )
        }
    }

    /// Apply a Hadamard transformation to `self` in place, swapping `0 ↔ +` and
    /// `1 ↔ -`.
    pub fn hadamard_mut(&mut self) {
        match *self {
            Self::Zero  => { *self = Self::Plus  },
            Self::One   => { *self = Self::Minus },
            Self::Plus  => { *self = Self::Zero  },
            Self::Minus => { *self = Self::One   },
        }
    }

    /// Apply a Hadamard transformation `self`, swapping `0 ↔ +` and `1 ↔ -`.
    pub fn hadamard(&self) -> Self {
        match self {
            Self::Zero  => Self::Plus,
            Self::One   => Self::Minus,
            Self::Plus  => Self::Zero,
            Self::Minus => Self::One,
        }
    }
}

impl StateOrd for State {
    fn state_cmp(&self, other: &Self) -> std::cmp::Ordering { self.cmp(other) }
}

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::Zero  => "0".fmt(f),
            Self::One   => "1".fmt(f),
            Self::Plus  => "+".fmt(f),
            Self::Minus => "-".fmt(f),
        }
    }
}

/// A specific basis for representing states.
///
/// The *Z* basis refers to [`State`]s `Zero` and `One`, while the *X* basis
/// refers to `Plus` and `Minus`. The *Z* basis is chosen by default for actual
/// computations.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Basis {
    Z,
    X,
}

impl Basis {
    // get the "up" (0/+) state for the basis
    pub(crate) fn up(&self) -> State {
        match self {
            Self::Z => State::Zero,
            Self::X => State::Plus,
        }
    }

    // get the "down" (1/-) state for the basis
    pub(crate) fn down(&self) -> State {
        match self {
            Self::Z => State::One,
            Self::X => State::Minus,
        }
    }
}

/// Describes a non-contiguously indexed collection of [`State`]s, to be used as
/// one "half" of a ketbra.
#[derive(Clone, Debug)]
pub struct States {
    pub(crate) offs: usize,
    pub(crate) data: VecDeque<Option<State>>,
    pub(crate) count: usize,
}

impl FromIterator<(usize, State)> for States {
    fn from_iter<I>(iter: I) -> Self
    where I: IntoIterator<Item = (usize, State)>
    {
        let mut data: VecDeque<Option<State>> = VecDeque::new();
        let mut k_min = usize::MAX;
        let mut k_max = usize::MIN;
        let mut count: usize = 0;
        for (k, s) in iter.into_iter() {
            if data.is_empty() {
                data.push_back(Some(s));
                k_min = k;
                k_max = k;
                count += 1;
            } else if (k_min..=k_max).contains(&k) {
                if data[k - k_min].is_none() { count += 1; }
                data[k - k_min] = Some(s);
            } else if k < k_min {
                (0..k_min - k - 1).for_each(|_| { data.push_front(None); });
                data.push_front(Some(s));
                k_min = k;
                count += 1;
            } else if k > k_max {
                (0..k - k_max - 1).for_each(|_| { data.push_back(None); });
                data.push_back(Some(s));
                k_max = k;
                count += 1;
            } else {
                unreachable!()
            }
        }
        Self { offs: k_min, data, count }
    }
}

impl IntoIterator for States {
    type Item = (usize, State);
    type IntoIter = StatesIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        StatesIntoIter {
            offs: self.offs,
            len: self.count,
            iter: self.data.into_iter().enumerate(),
        }
    }
}

impl Default for States {
    fn default() -> Self { Self::new() }
}

impl States {
    /// Create a new, empty `States`.
    pub fn new() -> Self {
        Self { offs: 0, data: VecDeque::new(), count: 0 }
    }

    /// Create a new, empty `States` with pre-allocated, contiguous indices over
    /// the range `a..=b`.
    pub fn with_capacity(a: usize, b: usize) -> Self {
        let kmin = a.min(b);
        let kmax = a.max(b);
        let mut data: VecDeque<Option<State>> =
            VecDeque::with_capacity(kmax - kmin + 1); // might be inefficient
        (kmin..=kmax).for_each(|_| { data.push_back(None); });
        data.make_contiguous();
        Self { offs: kmin, data, count: 0 }
    }

    /// Return the maximal `a` and minimal `b` such that all inhabited indices
    /// fall in the range `a..=b`.
    ///
    /// Returns `None` if no indices are inhabited.
    pub fn idx_bounds(&self) -> Option<(usize, usize)> {
        let kmin =
            self.data.iter().enumerate()
            .find_map(|(id, mb_s)| mb_s.map(|_| id + self.offs))?;
        let kmax =
            self.data.iter().enumerate()
            .rfind(|(_, mb_s)| mb_s.is_some())
            .map(|(id, _)| id + self.offs)?;
        Some((kmin, kmax))
    }

    /// Shift all indices to the right (larger values) by a uniform distance.
    ///
    /// *Panics if the size of the shift would cause any indices to overflow
    /// `usize`.*
    pub fn shift_indices_r(&mut self, sh: usize) {
        if let Some((_, kmax)) = self.idx_bounds() {
            if usize::MAX - kmax - 1 < sh {
                panic!("States::shift_indices_r: shift would cause index overflow");
            }
        }
        self.offs += sh;
    }

    /// Shift all indices to the left (smaller values) by a uniform distance.
    ///
    /// *Panics if the size of the shift would cause any indices to underflow
    /// `usize`.*
    pub fn shift_indices_l(&mut self, sh: usize) {
        if let Some((kmin, _)) = self.idx_bounds() {
            if kmin < sh {
                panic!("States::shift_indices_l: shift would cause index underflow");
            }
        }
        if self.offs < sh {
            (0..sh - self.offs).for_each(|_| { self.data.pop_front(); });
            self.offs = 0;
        } else {
            self.offs -= sh;
        }
    }

    /// Create a new `States` containing a single state.
    pub fn single(k: usize, state: State) -> Self {
        let mut data: VecDeque<Option<State>> = VecDeque::new();
        data.push_back(Some(state));
        Self { offs: k, data, count: 1 }
    }

    /// Get a reference to the state associated with the `k`-th index.
    pub fn get(&self, k: usize) -> Option<&State> {
        if k < self.offs || k > self.offs + self.data.len() - 1 {
            None
        } else {
            self.data[k - self.offs].as_ref()
        }
    }

    /// Get a mutable reference to the state associated with the `k`-th index.
    pub fn get_mut(&mut self, k: usize) -> Option<&mut State> {
        if k < self.offs || k > self.offs + self.data.len() - 1 {
            None
        } else {
            self.data[k - self.offs].as_mut()
        }
    }

    /// Set the state of the `k`-th index.
    ///
    /// Existing states are replaced and returned, if they exist.
    pub fn insert(&mut self, k: usize, state: State) -> Option<State> {
        if self.data.is_empty() {
            self.offs = k;
            self.data.push_back(Some(state));
            self.count += 1;
            None
        } else if k < self.offs {
            (0..self.offs - k - 1)
                .for_each(|_| { self.data.push_front(None); });
            self.data.push_front(Some(state));
            self.offs = k;
            self.count += 1;
            None
        } else if k > self.offs + self.data.len() - 1 {
            (0..k - self.data.len() - self.offs)
                .for_each(|_| { self.data.push_back(None); });
            self.data.push_back(Some(state));
            self.count += 1;
            None
        } else {
            let ret = self.data[k - self.offs].replace(state);
            if ret.is_none() { self.count += 1; }
            ret
        }
    }

    /// Remove and return the `k`-th state from the collection, if it exists.
    pub fn remove(&mut self, k: usize) -> Option<State> {
        if (self.offs..self.offs + self.data.len()).contains(&k) {
            let ret = self.data[k - self.offs].take();
            if ret.is_some() { self.count -= 1; }
            ret
        } else {
            None
        }
    }

    /// Swap the `a`-th and `b`-th states.
    ///
    /// If only one of either exist, the existing state is moved to the new
    /// position. If neither exist, no operation is performed.
    pub fn swap(&mut self, a: usize, b: usize) {
        if a == b || !(self.contains_key(a) || self.contains_key(b)) { return; }
        let range = self.offs..self.offs + self.data.len();
        let kmin = a.min(b);
        let kmax = a.max(b);
        if !range.contains(&kmin) {
            (0..self.offs - kmin)
                .for_each(|_| { self.data.push_front(None); });
            self.offs = kmin;
        }
        if !range.contains(&kmax) {
            (0..kmax - self.data.len() - self.offs + 1)
                .for_each(|_| { self.data.push_back(None); });
        }
        println!("{} {} {:?}", kmin, kmax, self.data);
        self.data.swap(a - self.offs, b - self.offs);
    }

    pub(crate) fn make_contiguous(&mut self) {
        self.data.make_contiguous();
    }

    /// Return the number of occupied indices.
    pub fn len(&self) -> usize { self.count }

    /// Return `true` if no indices are occupied.
    pub fn is_empty(&self) -> bool { self.count == 0 }

    /// Return an iterator over all occupied indices.
    ///
    /// The iterator item type is `(usize, &`[`State`]`)`.
    pub fn iter(&self) -> StatesIter<'_> {
        StatesIter {
            offs: self.offs,
            len: self.count,
            iter: self.data.iter().enumerate(),
        }
    }

    /// Return an iterator over mutable references to all occupied indices.
    ///
    /// The iterator item type is `(usize, &mut `[`State`]`)`.
    pub fn iter_mut(&mut self) -> StatesIterMut<'_> {
        StatesIterMut {
            offs: self.offs,
            len: self.count,
            iter: self.data.iter_mut().enumerate(),
        }
    }

    /// Return `true` if the `k`-th state exists in the collection.
    pub fn contains_key(&self, k: usize) -> bool {
        if (self.offs..self.offs + self.data.len()).contains(&k) {
            self.data[k - self.offs].is_some()
        } else {
            false
        }
    }
}

impl fmt::Display for States {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[")?;
        let n = self.len();
        for (k, (id, s)) in self.iter().enumerate() {
            write!(f, "{}: ", id)?;
            s.fmt(f)?;
            if k < n - 1 { write!(f, ", ")?; }
        }
        write!(f, "]")?;
        Ok(())
    }
}

impl PartialEq for States {
    fn eq(&self, other: &Self) -> bool {
        if self.count != other.count {
            false
        } else {
            self.iter().zip(other.iter())
                .all(|((id_l, s_l), (id_r, s_r))| id_l == id_r && s_l == s_r)
        }
    }
}

impl Eq for States { }

impl PartialOrd for States {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for States {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use itertools::EitherOrBoth::*;
        for mb_pair in self.iter().zip_longest(other.iter()) {
            match mb_pair {
                Right(_) => { return std::cmp::Ordering::Less; },
                Left(_) => { return std::cmp::Ordering::Greater; },
                Both((id_l, s_l), (id_r, s_r)) => {
                    if id_l < id_r {
                        return std::cmp::Ordering::Less;
                    } else if id_l > id_r {
                        return std::cmp::Ordering::Greater;
                    } else if s_l < s_r {
                        return std::cmp::Ordering::Less;
                    } else if s_l > s_r {
                        return std::cmp::Ordering::Greater;
                    } else {
                        continue;
                    }
                },
            }
        }
        std::cmp::Ordering::Equal
    }
}

impl StateOrd for States {
    fn state_cmp(&self, other: &Self) -> std::cmp::Ordering { self.cmp(other) }
}

/// Iterator over the occupied indices of a [`States`].
///
/// The iterator element type is `(usize, &`[`State`]`)`.
#[derive(Clone, Debug)]
pub struct StatesIter<'a> {
    offs: usize,
    len: usize,
    iter: std::iter::Enumerate<std::collections::vec_deque::Iter<'a, Option<State>>>,
}

impl<'a> Iterator for StatesIter<'a> {
    type Item = (usize, &'a State);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find_map(|(id, mb_s)| {
            mb_s.as_ref()
                .map(|s| {
                    self.len = self.len.saturating_sub(1);
                    (self.offs + id, s)
                })
        })
    }
}

impl<'a> DoubleEndedIterator for StatesIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(_, mb_s)| mb_s.is_some())
            .map(|(id, mb_s)| {
                self.len = self.len.saturating_sub(1);
                (self.offs + id, mb_s.as_ref().unwrap())
            })
    }
}

impl<'a> ExactSizeIterator for StatesIter<'a> {
    fn len(&self) -> usize { self.len }
}

impl<'a> std::iter::FusedIterator for StatesIter<'a> { }

/// Iterator over mutable references to the occupied indices of a [`States`].
///
/// The iterator element type is `(usize, &mut `[`State`]`)`.
#[derive(Debug)]
pub struct StatesIterMut<'a> {
    offs: usize,
    len: usize,
    iter: std::iter::Enumerate<std::collections::vec_deque::IterMut<'a, Option<State>>>,
}

impl<'a> Iterator for StatesIterMut<'a> {
    type Item = (usize, &'a mut State);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find_map(|(id, mb_s)| {
            mb_s.as_mut()
                .map(|s| {
                    self.len = self.len.saturating_sub(1);
                    (self.offs + id, s)
                })
        })
    }
}

impl<'a> DoubleEndedIterator for StatesIterMut<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(_, mb_s)| mb_s.is_some())
            .map(|(id, mb_s)| {
                self.len = self.len.saturating_sub(1);
                (self.offs + id, mb_s.as_mut().unwrap())
            })
    }
}

impl<'a> ExactSizeIterator for StatesIterMut<'a> {
    fn len(&self) -> usize { self.len }
}

impl<'a> std::iter::FusedIterator for StatesIterMut<'a> { }

/// Iterator over all occupied indices of a [`States`].
///
/// The iterator element type is `(usize, `[`State`]`)`.
#[derive(Clone, Debug)]
pub struct StatesIntoIter {
    offs: usize,
    len: usize,
    iter: std::iter::Enumerate<std::collections::vec_deque::IntoIter<Option<State>>>,
}

impl Iterator for StatesIntoIter {
    type Item = (usize, State);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.find_map(|(id, mb_s)| {
            mb_s.map(|s| {
                self.len = self.len.saturating_sub(1);
                (self.offs + id, s)
            })
        })
    }
}

impl DoubleEndedIterator for StatesIntoIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.rfind(|(_, mb_s)| mb_s.is_some())
            .map(|(id, mb_s)| {
                self.len = self.len.saturating_sub(1);
                (self.offs + id, mb_s.unwrap())
            })
    }
}

impl ExactSizeIterator for StatesIntoIter {
    fn len(&self) -> usize { self.len }
}

impl std::iter::FusedIterator for StatesIntoIter { }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decomp() {
        use Basis::*;
        use State::*;
        let rt2 = C64::from(std::f64::consts::FRAC_1_SQRT_2);
        assert_eq!(Zero.decomp(Z),  None                              );
        assert_eq!(Zero.decomp(X),  Some(((rt2, Plus), ( rt2, Minus))));
        assert_eq!(One.decomp(Z),   None                              );
        assert_eq!(One.decomp(X),   Some(((rt2, Plus), (-rt2, Minus))));
        assert_eq!(Plus.decomp(Z),  Some(((rt2, Zero), ( rt2, One)  )));
        assert_eq!(Plus.decomp(X),  None                              );
        assert_eq!(Minus.decomp(Z), Some(((rt2, Zero), (-rt2, One)  )));
        assert_eq!(Minus.decomp(X), None                              );
    }

    #[test]
    fn from_iter() {
        let states: States =
            [
                (2, State::Zero),
                (5, State::Plus),
                (3, State::One),
                (8, State::Minus),
            ]
            .into_iter()
            .collect();
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, 2);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 4);

        let states: States =
            [
                (1, State::Zero),
                (1, State::Plus),
                (1, State::Minus),
                (1, State::One),
            ]
            .into_iter()
            .collect();
        let data_expected: VecDeque<Option<State>> =
            vec![Some(State::One)].into();
        assert_eq!(states.offs, 1);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 1);
    }

    #[test]
    fn with_capacity() {
        let states = States::with_capacity(5, 7);
        let data_expected: VecDeque<Option<State>> =
            vec![None, None, None].into();
        assert_eq!(states.offs, 5);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 0);
    }

    #[test]
    fn insert() {
        let mut states = States::new();
        states.insert(2, State::Zero);
        states.insert(5, State::Plus);
        states.insert(3, State::One);
        states.insert(8, State::Minus);
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, 2);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 4);
    }

    fn build_simple() -> States {
        let mut states = States::new();
        states.insert(2, State::Zero);
        states.insert(1, State::Zero);
        states.insert(5, State::Plus);
        states.insert(3, State::One);
        states.insert(8, State::Minus);
        states
    }

    #[test]
    fn remove() {
        let mut states = build_simple();
        assert_eq!(states.remove(4), None);
        assert_eq!(states.remove(9), None);
        assert_eq!(states.remove(0), None);
        assert_eq!(states.remove(3), Some(State::One));
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Zero),
                None,
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, 1);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 4);
    }

    #[test]
    fn idx_bounds() {
        let mut states = build_simple();
        assert_eq!(states.idx_bounds(), Some((1, 8)));
        states.remove(2);
        assert_eq!(states.idx_bounds(), Some((1, 8)));
        states.remove(1);
        assert_eq!(states.idx_bounds(), Some((3, 8)));
        states.remove(8);
        assert_eq!(states.idx_bounds(), Some((3, 5)));
        states.remove(3);
        states.remove(5);
        assert_eq!(states.idx_bounds(), None);
    }

    #[test]
    fn shift() {
        let mut states = build_simple();
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, 1);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 5);

        states.shift_indices_l(1);
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, 0);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 5);

        states.shift_indices_r(5);
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, 5);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 5);

        states.insert(1, State::One);
        states.remove(1);
        let data_expected: VecDeque<Option<State>> =
            vec![
                None,
                None,
                None,
                None,
                Some(State::Zero),
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, 1);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 5);

        states.shift_indices_l(3);
        let data_expected: VecDeque<Option<State>> =
            vec![
                None,
                None,
                Some(State::Zero),
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, 0);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 5);

        states.shift_indices_l(2);
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, 0);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 5);

        states.shift_indices_r(usize::MAX - 8);
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.offs, usize::MAX - 8);
        assert_eq!(states.data, data_expected);
        assert_eq!(states.count, 5);

        assert!(states.get(usize::MAX - 1).is_some());
    }

    #[test]
    #[should_panic]
    fn bad_shift_l() {
        let mut states = build_simple();
        states.shift_indices_l(2);
    }

    #[test]
    #[should_panic]
    fn bad_shift_r() {
        let mut states = build_simple();
        states.shift_indices_r(usize::MAX - 6);
    }

    #[test]
    fn get() {
        let states = build_simple();
        assert_eq!(states.get(0), None);
        assert_eq!(states.get(1), Some(&State::Zero));
        assert_eq!(states.get(4), None);
        assert_eq!(states.get(10), None);
    }

    #[test]
    fn contains_key() {
        let states = build_simple();
        assert_eq!(states.contains_key(0), false);
        assert_eq!(states.contains_key(1), true);
        assert_eq!(states.contains_key(4), false);
        assert_eq!(states.contains_key(10), false);
    }

    #[test]
    fn eq() {
        let states0 = build_simple();
        let mut states1 = build_simple();
        states1.insert(10, State::Plus);
        states1.remove(10);
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                Some(State::Minus),
                None,
                None,
            ]
            .into();
        assert_eq!(states1.offs, 1);
        assert_eq!(states1.data, data_expected);
        assert_eq!(states1.count, 5);
        assert_eq!(states0, states1);
    }

    #[test]
    fn swap() {
        let mut states = build_simple();
        states.swap(8, 10);
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Zero),
                Some(State::One),
                None,
                Some(State::Plus),
                None,
                None,
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.data, data_expected);
        states.swap(2, 5);
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Plus),
                Some(State::One),
                None,
                Some(State::Zero),
                None,
                None,
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.data, data_expected);
        states.swap(0, 20);
        let data_expected: VecDeque<Option<State>> =
            vec![
                Some(State::Zero),
                Some(State::Plus),
                Some(State::One),
                None,
                Some(State::Zero),
                None,
                None,
                None,
                None,
                Some(State::Minus),
            ]
            .into();
        assert_eq!(states.data, data_expected);
    }

    #[test]
    fn len_is_empty() {
        let mut states = build_simple();
        assert_eq!(states.len(), 5);
        assert_eq!(states.is_empty(), false);
        states.remove(1);
        assert_eq!(states.len(), 4);
        assert_eq!(states.is_empty(), false);
        states.remove(0);
        assert_eq!(states.len(), 4);
        assert_eq!(states.is_empty(), false);
        states.remove(2);
        assert_eq!(states.len(), 3);
        assert_eq!(states.is_empty(), false);
        states.remove(3);
        assert_eq!(states.len(), 2);
        assert_eq!(states.is_empty(), false);
        states.remove(4);
        assert_eq!(states.len(), 2);
        assert_eq!(states.is_empty(), false);
        states.remove(5);
        states.remove(8);
        assert_eq!(states.len(), 0);
        assert_eq!(states.is_empty(), true);
    }

    #[test]
    fn iter() {
        let states = build_simple();
        let items: Vec<(usize, &State)> = states.iter().collect();
        let items_expected: Vec<(usize, &State)> =
            vec![
                (1, &State::Zero),
                (2, &State::Zero),
                (3, &State::One),
                (5, &State::Plus),
                (8, &State::Minus),
            ];
        assert_eq!(items, items_expected);
    }

    #[test]
    fn into_iter() {
        let states = build_simple();
        let items: Vec<(usize, State)> = states.into_iter().collect();
        let items_expected: Vec<(usize, State)> =
            vec![
                (1, State::Zero),
                (2, State::Zero),
                (3, State::One),
                (5, State::Plus),
                (8, State::Minus),
            ];
        assert_eq!(items, items_expected);
    }

}

