use nalgebra::Vector3;

#[derive(PartialEq)]
#[derive(Debug)]
pub struct Mission {
    pub kepler_elements: KeplerElements,
    pub ode_solver: OdeSolver,
    pub epoch: f64,
    pub coordinate_system: CoordinateSystem,
}

impl Default for Mission {
    fn default() -> Self {
        Mission {
            kepler_elements: KeplerElements {
                ..Default::default()
            },
            ode_solver: OdeSolver::RungeKutta4,
            epoch: 0.0,
            coordinate_system: CoordinateSystem::EarthCenteredInertial,
        }
    }
}

#[derive(Debug)]
#[derive(PartialEq)]
pub enum OdeSolver {
    RungeKutta4,
    DormandPrince5,
    DormandPrince853,
}

#[derive(Debug)]
#[derive(PartialEq)]
pub enum CoordinateSystem {
    EarthCenteredInertial,
    EarthCenteredEarthFixed,
}

///! **Kepler Elements**
///! Struct which defines the basic Kepler elements
#[derive(PartialEq)]
#[derive(Debug)]
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
            r * (self.longitude_of_ascending_node.cos() * (arg_of_lat).cos()
                - self.longitude_of_ascending_node.sin()
                    * (arg_of_lat).sin()
                    * self.inclination.cos()),
            r * (self.longitude_of_ascending_node.sin() * (arg_of_lat).cos()
                + self.longitude_of_ascending_node.sin()
                    * (arg_of_lat).sin()
                    * self.inclination.cos()),
            r * (self.argument_of_periapsis.cos() * self.inclination.sin()),
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
        let kepler_elements = KeplerElements {
            semi_major_axis: 7000.0,
            eccentricity: 0.1,
            inclination: 0.0,
            longitude_of_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            true_anomaly: 0.0,
        };

        let mu = 398600.4418e9;

        let (position, velocity) = kepler_elements.to_state_vector(mu);

        assert_eq!(position, Vector3::new(6299.999999999999, 0.0, 0.0));
        assert_eq!(velocity, Vector3::new(0.0, 263812.24864760914, 0.0));
    }
}