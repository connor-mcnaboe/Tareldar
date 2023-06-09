use std::fmt;
use std::fmt::Formatter;
use std::str::FromStr;

// TODO: Figure out a better way to load body data, this is a hacky quick solution.
pub fn get_body(body: &CentralBody) -> Body {
    match body {
        CentralBody::EARTH => Body { mu: 398600.4418E9 }
    }
}

pub struct Body {
    pub mu: f64,
}

#[derive(PartialEq, Debug)]
pub enum CentralBody {
    EARTH,
}

impl fmt::Display for CentralBody {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            CentralBody::EARTH => write!(f, "EARTH"),
        }
    }
}

impl FromStr for CentralBody {
    type Err = ();

    fn from_str(input: &str) -> Result<CentralBody, Self::Err> {
        match input {
            "EARTH" => Ok(CentralBody::EARTH),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod core_tests {
    use super::*;

    #[test]
    fn test_central_body_enum_supports_to_string() {
        assert_eq!(CentralBody::EARTH.to_string(), "EARTH");
    }

    #[test]
    fn test_central_body_enum_supports_from_str() {
        assert_eq!(CentralBody::from_str("EARTH").unwrap(), CentralBody::EARTH);
    }
}
