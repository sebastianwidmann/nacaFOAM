from numpy import sqrt, power, log10

# Define thermodynamic properties of air at ICAO standard atmosphere
T0 = 288.15  # [K] Total temperature
p0 = 101325  # [Pa] Total pressure
gamma = 1.4  # [-] Ratio of specific heats
R = 287.058  # [J/(kg*K)] Specific gas constant for dry air

mu = 1.789e-5  # [] Dynamic viscosity
cp = 1005  # [J/(kg*K)
lamb = 25.5e-3  # [W/(m*K) thermal conductivity


def calculateStaticTemperature(M):
    return T0 / (1 + 0.5 * (gamma - 1) * M ** 2)


def calculateStaticPressure(M):
    return p0 * (1 + 0.5 * (gamma - 1) * M ** 2) ** (-gamma / (gamma - 1))


def calculateStaticDensity(p, T):
    return p / (R * T)


def calculateSpeedofSound(T):
    return sqrt(gamma * R * T)


def calculateDynamicViscosity(T):
    return 1.716e-5 * (T / 273.15) ** (1.5) * (273.15 + 110.4) / (T + 110.4)


def calculatePrandtlNumber():
    return mu * cp / lamb

def calculateReynoldsNumber(rho, u, L, mu):
    return rho * u * L / mu

def calculateFirstLayerThickness(mach, yPlus):
    """
    Returns
    -------
    Minimum cell height in boundary layer
    """

    p = calculateStaticPressure(mach)
    T = calculateStaticTemperature(mach)
    rho = calculateStaticDensity(p, T)
    u = mach * calculateSpeedofSound(T)
    Re = rho * u / mu  # [-] Freestream Reynolds number
    cf = power(2 * log10(Re) - 0.65, -2.3)  # [-] Skin friction coefficient based on Schlichting
    Tau_w = 0.5 * cf * rho * u ** 2  # [Pa] Wall shear stress
    u_star = sqrt(Tau_w / rho)  # [m*s^-1] Friction velocity
    return yPlus * mu / (rho * u_star)
