from heat_transfer import piping as pipe
import logging

ureg = pipe.ureg
Q_ = ureg.Quantity

logger = logging.getLogger(__name__)


def min_thick(P, T, Pipe, S=16700*ureg.psi, E=0.8, W=1, c=0*ureg.inch):
    """Calculate min required thickness for Straight Pipe Under Internal Pressure as per B31.3 304.1.2
    :param P: Design pressure
    :param T: Design temperature
    :param Pipe: Pipe object
    :param S: Allowable stress from table A-1
    :param E: Quality factor from table A-1A, A-1B
    :param W: Weld joint strength reduction factor from table 302.3.5
    :return: return description
    :rtype: the return type description
    """
    D = Pipe.OD
    Y = Y_coef(Pipe.ID, Pipe.OD, T, c)
    t = P*D/(2*(S*E*W)+P*Y)
    return t+c

def Y_coef(d, D, T, c=0*ureg.inch):
    if (D-d)/2<D/6:
        if T<Q_(482, ureg.degC): #TODO add interpolation of the table 304.1.1
            Y = 0.4
        else:
            logger.error('Temperatures above 482 degC are not implemented yet: {}'.format(T))
    else:
        Y = (d+2*c)/(D+d+2*c)
    return Y




Test_pipe = pipe.Pipe(1, SCH=10)
print(Test_pipe.OD)
P = ureg("100psi")
T = ureg("100K")
print(min_thick(P, T, Test_pipe))
