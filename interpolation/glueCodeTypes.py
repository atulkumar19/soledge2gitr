from enum import Enum, IntEnum
import collections

class ALInterfaceMode(IntEnum):
    FGS = 0
    ACTIVELEARNER=2
    FAKE = 3
    DEFAULT = 4
    FASTFGS = 5
    ANALYTIC = 6
    KILL = 9

class SolverCode(Enum):
    BGK = 0
    LBMZEROD = 1
    BGKMASSES = 2

class DatabaseMode(IntEnum):
    SQLITE = 0

# BGKInputs

BGKInputs = collections.namedtuple('BGKInputs', 'R Z')
BGKOutputs = collections.namedtuple('BGKOutputs', 'BPHI BR BZ TE NE VE TI NI VI')
