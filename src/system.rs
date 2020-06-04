use std::fmt::{Display, Formatter, Result as FmtResult};
use petgraph::{graph::{NodeIndex, DiGraph}, algo::tarjan_scc, dot::Dot};

type Ix = u32;

/// A variable in a `System`.
pub type Var = NodeIndex<Ix>;

/// A relation for the system of equations.
#[derive(Debug, Clone, Copy, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum Relation {
    /// Used for equality as well.
    LessEqual,
    Less
}

/// A system of equations.
#[derive(Debug, Clone)]
pub struct System(DiGraph<(), Relation, Ix>);

impl Default for System {
    fn default() -> Self {
        System::new()
    }
}

impl Display for System {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "{:?}", Dot::new(&self.0))
    }
}

impl System {
    /// Create a system with no ordering constraints or variables.
    pub fn new() -> Self {
        System(DiGraph::with_capacity(100, 100))
    }

    /// Creates a new variable for the constraint system.
    pub fn new_variable(&mut self) -> Var {
        self.0.add_node(())
    }
    
    /// Returns whether the system is satisfiable.
    ///
    /// See Theorem 1 in [*SMT Solving for the Theory of Ordering Constraints*] by
    /// Cunjing Ge, Feifei Ma, Jeff Huang, and Jian Zhang to see how this works.
    ///
    /// [*SMT Solving for the Theory of Ordering Constraints*]: https://parasol.tamu.edu/groups/huangroup/academic/coco.pdf
    pub fn is_satisfiable(&self) -> bool {
        for scc in tarjan_scc(&self.0) {
            for v1 in scc.iter() {
                for v2 in self.0.neighbors(*v1) {
                    if scc.contains(&v2) {
                        for e in self.0.edges_connecting(*v1, v2) {
                            if *e.weight() == Relation::Less {
                                return false
                            }
                        }
                    }
                }
            }
        }

        true
    }

    /// Inserts a new constraint to the system.
    #[inline]
    pub fn push(&mut self, l: Var, rel: Relation, r: Var) {
        self.0.add_edge(l, r, rel);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Should be satisfiable.
    #[test]
    fn big_system_satisfiable() {
        let mut sys = System::new();
        let mut vars = Vec::new();

        for _ in 0..5000 {
            let var = sys.new_variable();
            vars.push(var)
        }

        let mut le = true;

        for lr in vars.windows(2) {
            if let [l, r] = lr {
                if le {
                    sys.push(*l, Relation::LessEqual, *r)
                } else {
                    sys.push(*l, Relation::Less, *r)
                }
                
                le = !le;
            } else {
                panic!()
            }
        }

        assert!(sys.is_satisfiable())
    }

    /// Make sure System can't solve unsolvable systems.
    #[test]
    fn small_system_fail() {
        let mut sys = System::new();
        let var1 = sys.new_variable();
        let var2 = sys.new_variable();

        sys.push(var1, Relation::Less, var2);
        sys.push(var2, Relation::Less, var1);

        assert!(!sys.is_satisfiable());
    }

    /// Much bigger unsatisfiable system.
    #[test]
    fn big_system_fail() {
        let mut sys = System::new();
        let mut vars = Vec::with_capacity(100);

        for _ in 0..5000 {
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

        assert!(!sys.is_satisfiable())
    }
}
