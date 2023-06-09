use crate::bodies::get_body;
use crate::mission::Mission;
use crate::orbit::Orbit;
use ode_solvers::dopri5::*;
use ode_solvers::*;
use std::fmt;
use std::fmt::Formatter;
use std::str::FromStr;

type State = Vector6<f64>;
type Time = f64;

#[derive(Debug, PartialEq)]
pub enum OdeSolver {
    RungeKutta4,
    DormandPrince5,
    DormandPrince853,
}

impl fmt::Display for OdeSolver {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            OdeSolver::RungeKutta4 => write!(f, "RungeKutta4"),
            OdeSolver::DormandPrince5 => write!(f, "DormandPrince5"),
            OdeSolver::DormandPrince853 => write!(f, "DormandPrince853"),
        }
    }
}

impl FromStr for OdeSolver {
    type Err = ();

    fn from_str(input: &str) -> Result<OdeSolver, Self::Err> {
        match input {
            "RungeKutta4" => Ok(OdeSolver::RungeKutta4),
            "DormandPrince5" => Ok(OdeSolver::DormandPrince5),
            "DormandPrince853" => Ok(OdeSolver::DormandPrince853),
            _ => Err(()),
        }
    }
}

impl System<State> for Orbit {
    /**
    Kepler Orbit Equations of motion.
    # Arguments
       * '_t' - The moment in time corresponding to a specific state.
       * 'y' - The state vector
       * 'dy' -  The change in the state vector
     */
    fn system(&self, _t: Time, y: &State, dy: &mut State) {
        let denominator: f64 = (y[0].powf(2.0) + y[1].powf(2.0) + y[2].powf(2.0)).powf(3.0 / 2.0);
        let body = get_body(&self.central_body);
        dy[0] = y[3];
        dy[1] = y[4];
        dy[2] = y[5];
        dy[3] = -body.mu * y[0] / denominator;
        dy[4] = -body.mu * y[1] / denominator;
        dy[5] = -body.mu * y[2] / denominator;
    }
}

/**
Propagate a state vector for a given time of flight.
 */
pub fn propagate(mission: Mission) -> Vec<Vector6<f64>> {
    let body = get_body(&mission.orbit.central_body);
    let (position, velocity) = mission.orbit.kepler_elements.to_state_vector(body.mu);
    let system = mission.orbit;

    let rtol: f64 = 1e-6;
    let atol: f64 = 1e-8;

    let mut stepper = Dopri5::new(
        system,
        mission.epoch,
        mission.duration,
        10.0,
        State::new(
            position[0],
            position[1],
            position[2],
            velocity[0],
            velocity[1],
            velocity[2],
        ),
        rtol,
        atol,
    );
    stepper
        .integrate()
        .expect("ERROR: Unable to integrate provided parameters.");

    let y_out = stepper.y_out();
    y_out.to_vec()
}

#[cfg(test)]
mod propagator_tests {
    use super::*;
    use crate::bodies::CentralBody;
    use crate::orbit::{CoordinateSystem, KeplerElements};
    use std::f64::consts::PI;

    #[test]
    fn test_should_integrate() {
        let eps_pos = 200.0;
        let eps_vel = 10.0;
        let mission = Mission {
            orbit: Orbit {
                kepler_elements: KeplerElements {
                    semi_major_axis: 6.791301224674748E+06,
                    eccentricity: 8.510618198049622E-04,
                    inclination: 4.949314343620572E+01 * PI / 180.0,
                    longitude_of_ascending_node: 9.440099680297747E+01 * PI / 180.0,
                    argument_of_periapsis: 8.122131421322101E+01 * PI / 180.0,
                    true_anomaly: 3.244321752988205E+02 * PI / 180.0,
                },
                central_body: CentralBody::EARTH,
                coordinate_system: CoordinateSystem::EarthCenteredInertial,
                ode_solver: OdeSolver::DormandPrince5,
            },
            epoch: 0.0,
            duration: 60.0*60.0,
        };
        let result = propagate(mission);
        let final_value = result.last().unwrap();

        // TODO: These are inaccurate values that will need to be updated as the model is improved.
        let expected_out_state = [
            4278239.0,
            1324790.0,
            -5111879.0,
            -1305.0,
            7494.0,
            851.0,
        ];

        // Position Vectors
        assert_relatively_eq(final_value[0], expected_out_state[0], eps_pos);
        assert_relatively_eq(final_value[1], expected_out_state[1], eps_pos);
        assert_relatively_eq(final_value[2], expected_out_state[2], eps_pos);

        // Velocity Vectors
        assert_relatively_eq(final_value[3], expected_out_state[3], eps_vel);
        assert_relatively_eq(final_value[4], expected_out_state[4], eps_vel);
        assert_relatively_eq(final_value[5], expected_out_state[5], eps_vel);
    }

    fn assert_relatively_eq(num_one: f64, num_two: f64, epsilon: f64) {
        let diff = (num_two - num_one).abs();
        assert!(diff <= epsilon);
    }

    #[test]
    fn test_ode_solver_enum_supports_to_string() {
        assert_eq!(OdeSolver::RungeKutta4.to_string(), "RungeKutta4");
        assert_eq!(OdeSolver::DormandPrince5.to_string(), "DormandPrince5");
        assert_eq!(OdeSolver::DormandPrince853.to_string(), "DormandPrince853")
    }

    #[test]
    fn test_ode_solver_enum_supports_from_str() {
        assert_eq!(
            OdeSolver::from_str("RungeKutta4").unwrap(),
            OdeSolver::RungeKutta4
        );
        assert_eq!(
            OdeSolver::from_str("DormandPrince5").unwrap(),
            OdeSolver::DormandPrince5
        );
        assert_eq!(
            OdeSolver::from_str("DormandPrince853").unwrap(),
            OdeSolver::DormandPrince853
        );
    }
}
