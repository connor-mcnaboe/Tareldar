use crate::orbit::Orbit;

#[derive(PartialEq, Debug)]
pub struct Mission {
    pub orbit: Orbit,
    pub epoch: f64,
    pub duration: f64,
}

impl Default for Mission {
    fn default() -> Self {
        Mission {
            orbit: Orbit::default(),
            epoch: 0.0,
            duration: 0.0,
        }
    }
}

#[cfg(test)]
mod mission_tests {
    use super::*;
    use crate::bodies::CentralBody;
    use crate::orbit::{CoordinateSystem, KeplerElements};
    use crate::propagator::OdeSolver;

    #[test]
    fn test_mission_sets_default() {
        let expected_elements = Mission {
            orbit: Orbit {
                kepler_elements: KeplerElements {
                    semi_major_axis: 1000.0,
                    eccentricity: 0.0,
                    inclination: 0.0,
                    longitude_of_ascending_node: 0.0,
                    argument_of_periapsis: 0.0,
                    true_anomaly: 0.0,
                },
                central_body: CentralBody::EARTH,
                coordinate_system: CoordinateSystem::EarthCenteredInertial,
                ode_solver: OdeSolver::RungeKutta4,
            },
            epoch: 0.0,
            duration: 0.0,
        };
        let actual_elements = Mission::default();
        assert_eq!(actual_elements, expected_elements)
    }
}
