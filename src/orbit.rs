use crate::bodies::CentralBody;
use crate::propagator::OdeSolver;
use nalgebra::Vector3;
use std::fmt;
use std::fmt::Formatter;
use std::str::FromStr;

#[derive(PartialEq, Debug)]
pub struct Orbit {
    pub kepler_elements: KeplerElements,
    pub central_body: CentralBody,
    pub coordinate_system: CoordinateSystem,
    pub ode_solver: OdeSolver,
}

impl Default for Orbit {
    fn default() -> Self {
        Orbit {
            kepler_elements: KeplerElements::default(),
            central_body: CentralBody::EARTH,
            coordinate_system: CoordinateSystem::EarthCenteredInertial,
            ode_solver: OdeSolver::RungeKutta4,
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum CoordinateSystem {
    EarthCenteredInertial,
    EarthCenteredEarthFixed,
}

impl fmt::Display for CoordinateSystem {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            CoordinateSystem::EarthCenteredInertial => write!(f, "EarthCenteredInertial"),
            CoordinateSystem::EarthCenteredEarthFixed => write!(f, "EarthCenteredEarthFixed"),
        }
    }
}

impl FromStr for CoordinateSystem {
    type Err = ();

    fn from_str(input: &str) -> Result<CoordinateSystem, Self::Err> {
        match input {
            "EarthCenteredInertial" => Ok(CoordinateSystem::EarthCenteredInertial),
            "EarthCenteredEarthFixed" => Ok(CoordinateSystem::EarthCenteredEarthFixed),
            _ => Err(()),
        }
    }
}

///! **Kepler Elements**
///! Struct which defines the basic Kepler elements
#[derive(PartialEq, Debug)]
pub struct KeplerElements {
    pub semi_major_axis: f64,             // in meters
    pub eccentricity: f64,                // dimensionless
    pub inclination: f64,                 // in radians
    pub longitude_of_ascending_node: f64, // in radians
    pub argument_of_periapsis: f64,       // in radians
    pub true_anomaly: f64,                // in radians
}

impl Default for KeplerElements {
    fn default() -> KeplerElements {
        KeplerElements {
            semi_major_axis: 1000.0,
            eccentricity: 0.0,
            inclination: 0.0,
            longitude_of_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            true_anomaly: 0.0,
        }
    }
}

impl KeplerElements {
    ///! Convert KeplerElements to position and velocity vectors
    pub fn to_state_vector(&self, mu: f64) -> (Vector3<f64>, Vector3<f64>) {
        //Calculate the Semi-latus rectum for the orbit:
        let semi_latus_rectum = self.semi_major_axis * (1.0 - self.eccentricity.powi(2));

        // Calculate the position and velocity in the orbital plane:
        let r = semi_latus_rectum / (1.0 + self.eccentricity * self.true_anomaly.cos());

        // Calculate the position and velocity in the Cartesian coordinate system
        let arg_of_lat = self.argument_of_periapsis + self.true_anomaly;
        let position = Vector3::new(
            r * ((arg_of_lat).cos() * self.longitude_of_ascending_node.cos()
                - (arg_of_lat).sin()
                    * self.inclination.cos()
                    * self.longitude_of_ascending_node.sin()),
            r * ((arg_of_lat).cos() * self.longitude_of_ascending_node.sin()
                + (arg_of_lat).sin()
                    * self.inclination.cos()
                    * self.longitude_of_ascending_node.cos()),
            r * ((arg_of_lat).sin() * self.inclination.sin()),
        );

        // Calculate the specific angular momentum vector
        let h = (mu * self.semi_major_axis * (1.0 - self.eccentricity.powi(2))).sqrt();

        let velocity = Vector3::new(
            ((position[0] * h * self.eccentricity) / (r * semi_latus_rectum))
                * self.true_anomaly.sin()
                - (h / r)
                    * (self.longitude_of_ascending_node.cos() * arg_of_lat.sin()
                        + self.longitude_of_ascending_node.sin()
                            * arg_of_lat.cos()
                            * self.inclination.cos()),
            ((position[1] * h * self.eccentricity) / (r * semi_latus_rectum))
                * self.true_anomaly.sin()
                - (h / r)
                    * (self.longitude_of_ascending_node.sin() * arg_of_lat.sin()
                        - self.longitude_of_ascending_node.cos()
                            * arg_of_lat.cos()
                            * self.inclination.cos()),
            ((position[1] * h * self.eccentricity) / (r * semi_latus_rectum))
                * self.true_anomaly.sin()
                + (h / r) * self.inclination.sin() * arg_of_lat.cos(),
        );

        (position, velocity)
    }
}

#[cfg(test)]
mod core_tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_orbit_sets_default() {
        let expected_elements = Orbit {
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
        };
        let actual_elements = Orbit::default();
        assert_eq!(actual_elements, expected_elements)
    }

    #[test]
    fn test_kepler_elements_sets_default() {
        let expected_elements = KeplerElements {
            semi_major_axis: 1000.0,
            eccentricity: 0.0,
            inclination: 0.0,
            longitude_of_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let actual_elements = KeplerElements {
            ..Default::default()
        };
        assert_eq!(actual_elements, expected_elements)
    }

    #[test]
    fn test_kepler_elements_to_state_vector() {
        let eps_pos = 0.5;
        let eps_vel = 0.5;
        let kepler_elements = KeplerElements {
            semi_major_axis: 6.791301224674748E+06,
            eccentricity: 8.510618198049622E-04,
            inclination: 4.949314343620572E+01 * PI / 180.0,
            longitude_of_ascending_node: 9.440099680297747E+01 * PI / 180.0,
            argument_of_periapsis: 8.122131421322101E+01 * PI / 180.0,
            true_anomaly: 3.244321752988205E+02 * PI / 180.0,
        };

        let expected_position = Vector3::new(
            -3.507115480698001E+06,
            4.487914768092333E+06,
            3.690078020656149E+06
        );

        let expected_velocity =  Vector3::new(
            -3.047824659988581E+03,
            -5.735901699645484E+03,
            4.072387230323159E+03
        );

        let mu = 398600.4418e9;

        let (position, velocity) = kepler_elements.to_state_vector(mu);

        // Position Vectors
        assert_relatively_eq(position[0], expected_position[0], eps_pos);
        assert_relatively_eq(position[1], expected_position[1], eps_pos);
        assert_relatively_eq(position[2], expected_position[2], eps_pos);

        // Velocity Vectors
        assert_relatively_eq(velocity[0], expected_velocity[0], eps_vel);
        assert_relatively_eq(velocity[1], expected_velocity[1], eps_vel);
        assert_relatively_eq(velocity[2], expected_velocity[2], eps_vel);
    }

    fn assert_relatively_eq(num_one: f64, num_two: f64, epsilon: f64) {
        let diff = (num_two - num_one).abs();
        assert!(diff <= epsilon);
    }

    #[test]
    fn test_coordinate_system_enum_supports_to_string() {
        assert_eq!(
            CoordinateSystem::EarthCenteredInertial.to_string(),
            "EarthCenteredInertial"
        );
        assert_eq!(
            CoordinateSystem::EarthCenteredEarthFixed.to_string(),
            "EarthCenteredEarthFixed"
        );
    }

    #[test]
    fn test_coordinate_system_enum_supports_from_str() {
        assert_eq!(
            CoordinateSystem::from_str("EarthCenteredEarthFixed").unwrap(),
            CoordinateSystem::EarthCenteredEarthFixed
        );
        assert_eq!(
            CoordinateSystem::from_str("EarthCenteredInertial").unwrap(),
            CoordinateSystem::EarthCenteredInertial
        );
    }
}
