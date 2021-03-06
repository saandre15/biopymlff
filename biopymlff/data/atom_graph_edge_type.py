from enum import IntEnum

class AtomGraphEdgeType(IntEnum):
    SINGLE=1
    DOUBLE=2
    TRIPLE=3
    AMIDE=4
    AROMATIC=5
    DUMMY=6
    UNKNOWN=7
    NOT_CONNECTED=8