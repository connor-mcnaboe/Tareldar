use nalgebra::Vector3;

/// **Kepler Elements**
/// Struct which defines the basic Kepler elements
///
pub struct KeplerElements {
    pub semi_major_axis: f64,             // in meters
    pub eccentricity: f64,                // dimensionless
    pub inclination: f64,                 // in radians
    pub longitude_of_ascending_node: f64, // in radians
    pub argument_of_periapsis: f64,       // in radians
    pub mean_anomaly: f64,                // in radians
}

impl KeplerElements {
    pub fn to_state_vector(&self, mu: f64) -> (Vector3<f64>, Vector3<f64>) {
        // Compute the mean motion (n):
        let n = (mu / self.semi_major_axis.powi(3)).sqrt()

        let true_anomaly =

        // Calculate the specific angular momentum vector
        let h = (mu * self.semi_major_axis * (1.0 - self.eccentricity.powi(2))).sqrt();

        // Calculate the position and velocity in the orbital plane
        let r = self.semi_major_axis * (1.0 - self.eccentricity.powi(2))
            / (1.0 + self.eccentricity * self.mean_anomaly.cos());
        let v = (mu / h).sqrt() * self.eccentricity.sin();

        // Calculate the position and velocity in the Cartesian coordinate system
        let p = Vector3::new(
            r * (self.argument_of_periapsis + self.mean_anomaly).cos(),
            r * (self.argument_of_periapsis + self.mean_anomaly).sin(),
            0.0,
        );
        let v = Vector3::new(
            -v * (self.argument_of_periapsis + self.mean_anomaly).sin(),
            v * (self.argument_of_periapsis + self.mean_anomaly).cos(),
            0.0,
        ) + h / r
            * Vector3::new(
                -((self.argument_of_periapsis).cos() * (self.inclination).sin()),
                (self.argument_of_periapsis).sin() * (self.inclination).sin(),
                (self.inclination).cos(),
            );

        // Return the position and velocity vectors
        (p, v)
    }
}

#[cfg(test)]
mod core_tests {
    use super::*;

    #[test]
    fn test_kepler_elements_to_state_vector() {
        let kepler_elements = KeplerElements {
            semi_major_axis: 7000.0,
            eccentricity: 0.1,
            inclination: 0.0,
            longitude_of_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            mean_anomaly: 0.0,
        };

        let mu = 398600.4418e9;

        let (position, velocity) = kepler_elements.to_state_vector(mu);

        assert_eq!(position, Vector3::new(7000.0, 0.0, 0.0));
        assert_eq!(velocity, Vector3::new(0.0, 7441.958657728055, 0.0));
    }
}
