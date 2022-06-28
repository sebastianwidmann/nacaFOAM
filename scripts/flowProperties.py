from math import sqrt

# Define thermodynamic properties of air at ICAO standard atmosphere
T0 = 288.15  # [K] Total temperature
p0 = 101325  # [Pa] Total pressure
gamma = 1.4  # [-] Ratio of specific heats
R = 287.058  # [J/(kg*K)] Specific gas constant for dry air

mu = 1.789e-5  # [] Dynamic viscosity
cp = 1005  # [J/(kg*K)
lamb = 25.5 - 3  # [W/(m*K) thermal conductivity


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
