use nalgebra::{DMatrix, DVector};
use approx::relative_eq;
use super::R;

/// A variable in a `System`.
#[derive(Debug, Clone, Copy, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub struct Var(usize);

/// A relation for the system of equations.
#[derive(Debug, Clone, Copy, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum Relation {
    /// Used for equality as well.
    LessEqual,
    Less
}

/// A system of equations.
#[derive(Debug, Clone, PartialEq)]
pub struct System {
    constraints: DMatrix<R>,
    constants: DVector<R>
}

impl Default for System {
    fn default() -> Self {
        System::new()
    }
}

impl System {
    /// Create a system with no constraints or variables.
    pub fn new() -> Self {
        System {
            constraints: DMatrix::from_element(0, 0, 0.0),
            constants: DVector::from_element(0, 0.0)
        }
    }

    /// Creates a new variable for the linear system.
    pub fn new_variable(&mut self) -> Var {
        let i = self.constraints.ncols() + 1;
        self.constraints.resize_horizontally_mut(i, 0.0);
        Var(i)
    }
    
    /// Returns whether the system is at all satisfiable,
    pub fn is_satisfiable(&self, eps: R) -> bool {
        if let Ok(pinv) = self.constraints.clone().pseudo_inverse(eps) {
            relative_eq!(
                self.constraints.clone() * pinv * self.constants.clone(),
                self.constants,
                epsilon = eps
            )
        } else {
            false
        }
    }

    /// Removes redundant constraints from `self`.
    ///
    /// A constraint is redundant when there is the exact same constraint in the
    /// system modulo the constant. The constraint with the larger constant is kept.
    pub fn simplify(&mut self) {
        fn aux(sys: &mut System) -> bool {
            let rows = sys.constraints.nrows();

            for i in 0..rows {
                for j in i + 1..rows {
                    if sys.constraints.row(i) == sys.constraints.row(j) {
                        let i_c = sys.constants.get(i)
                            .expect("constraints and constants dont have same len")
                            .clone();
                        let j_c = sys.constants.get(i)
                            .expect("constraints and constants dont have same len")
                            .clone();
                        
                        let row_remove = if i_c < j_c {
                            i
                        } else {
                            j
                        };
                        
                        let mut temp_constraints = DMatrix::from_element(0, 0, 0.0);
                        let mut temp_constants = DVector::from_element(0, 0.0);

                        std::mem::swap(&mut temp_constraints, &mut sys.constraints);
                        std::mem::swap(&mut temp_constants, &mut sys.constants);
                        
                        std::mem::replace(
                            &mut sys.constraints,
                            temp_constraints.remove_row(row_remove)
                        );
                        std::mem::replace(
                            &mut sys.constants,
                            temp_constants.remove_row(row_remove)
                        );

                        return true;
                    }
                }
            }

            false
        }
        
        while aux(self) {}
    }

    /// Inserts a new constraint to the system.
    pub fn push(&mut self, l: Var, rel: Relation, r: Var) {
        let c = if rel == Relation::Less {
            nalgebra::one()
        } else {
            nalgebra::zero()
        };

        self.constants.extend(std::iter::once(c));
        let constraint_col = self.constraints.nrows() + 1;
        self.constraints.resize_vertically_mut(constraint_col, nalgebra::zero());
        std::mem::replace(
            self.constraints
                .get_mut((constraint_col - 1, l.0 - 1))
                .expect("bad left variable"),
            -1.0
        );
        std::mem::replace(
            self.constraints
                .get_mut((constraint_col - 1, r.0 - 1))
                .expect("bad right variable"),
            1.0
        );
    }

    /// Inserts new constraints into the system.
    ///
    /// The constraints are of the form `le <= r` and `lt < r` for every respective
    /// entry in `l_le` and `l_lt`.
    pub fn push_max(&mut self, l_le: Vec<Var>, l_lt: Vec<Var>, r: Var) {
        for l in l_le.into_iter() {
            self.push(l, Relation::LessEqual, r)
        }

        for l in l_lt.into_iter() {
            self.push(l, Relation::Less, r)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Assert the sizes of a fresh System's matrices.
    /// 
    /// * `System.constraints` should have 0 cols, 0 rows.
    /// * `System.constants` should have 1 col, 0 rows.
    #[test]
    fn empty_size() {
        let sys = System::new();

        assert_eq!(sys.constants.ncols(), 1);
        assert_eq!(sys.constants.nrows(), 0);
        assert_eq!(sys.constraints.ncols(), 0);
        assert_eq!(sys.constraints.nrows(), 0);
    }
    
    /// Assert the sizes of a System's matrices when making new vars.
    ///
    /// After insertion of the `i`th variable,
    /// * `System.constraints` should have `i` cols, 0 rows.
    /// * `System.constants` should have 1 col, 0 rows.
    #[test]
    fn new_var_size() {
        let mut sys = System::new();

        for i in 0..100 {
            assert_eq!(sys.constants.ncols(), 1);
            assert_eq!(sys.constants.nrows(), 0);
            assert_eq!(sys.constraints.ncols(), i);
            assert_eq!(sys.constraints.nrows(), 0);
            sys.new_variable();
        }
    }

    /// Assert the sizes of a System's matrices when making new vars and constraints.
    ///
    /// After insertion of the `c`th constraint and with `v` variables,
    /// * `System.constraints` should have `v` cols, `c` rows.
    /// * `System.constants` should have 1 col, `c` rows.
    #[test]
    fn new_var_constraint_size() {
        let mut sys = System::new();
        let mut vars = Vec::with_capacity(100);

        for _ in 0..100 {
            let var = sys.new_variable();
            vars.push(var)
        }

        let mut le = true;
        let mut constraints = 0;

        for lr in vars.windows(2) {
            if let [l, r] = lr {
                if le {
                    sys.push(*l, Relation::LessEqual, *r)
                } else {
                    sys.push(*l, Relation::Less, *r)
                }
                
                le = !le;
                constraints += 1;

                assert_eq!(sys.constants.ncols(), 1);
                assert_eq!(sys.constants.nrows(), constraints);
                assert_eq!(sys.constraints.ncols(), 100);
                assert_eq!(sys.constraints.nrows(), constraints);
            } else {
                panic!()
            }
        }
    }

    /// Assert the sizes of a System's matrices when making new vars and constraints.
    ///
    /// After insertion of the `c`th constraint and with `v` variables,
    /// * `System.constraints` should have `v` cols, `c` rows.
    /// * `System.constants` should have 1 col, `c` rows.
    ///
    /// Should also be satisfiable.
    #[test]
    fn big_system_satisfiable() {
        let mut sys = System::new();
        let mut vars = Vec::new();

        for _ in 0..50 {
            let var = sys.new_variable();
            vars.push(var)
        }

        let mut le = true;
        let mut constraints = 0;

        for lr in vars.windows(2) {
            if let [l, r] = lr {
                if le {
                    sys.push(*l, Relation::LessEqual, *r)
                } else {
                    sys.push(*l, Relation::Less, *r)
                }
                
                le = !le;
                constraints += 1;

                assert_eq!(sys.constants.ncols(), 1);
                assert_eq!(sys.constants.nrows(), constraints);
                assert_eq!(sys.constraints.ncols(), 50);
                assert_eq!(sys.constraints.nrows(), constraints);
            } else {
                panic!()
            }
        }

        assert!(sys.is_satisfiable(1e-7))
    }

    /// Make sure System can't solve unsolvable systems.
    #[test]
    fn small_system_fail() {
        let mut sys = System::new();
        let var1 = sys.new_variable();
        let var2 = sys.new_variable();

        sys.push(var1, Relation::Less, var2);
        sys.push(var2, Relation::Less, var1);

        assert!(!sys.is_satisfiable(1e-7));
    }

    /// Much bigger unsatisfiable system.
    #[test]
    fn big_system_fail() {
        let mut sys = System::new();
        let mut vars = Vec::with_capacity(100);

        for _ in 0..50 {
            let var = sys.new_variable();
            vars.push(var)
        }

        for lr in vars.windows(2) {
            if let [l, r] = lr {
                sys.push(*l, Relation::Less, *r);
            } else {
                panic!()
            }
        }

        sys.push(*vars.last().unwrap(), Relation::LessEqual, *vars.first().unwrap());

        assert!(!sys.is_satisfiable(1e-7))
    }
}
