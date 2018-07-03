# Double Crystal Ball function multiplied by a constant for normalization
# p[0]-mean, p[1]-width, p[2]-alpha1, p[3]-n1, p[4]-alpha2, p[5]-n2, p[6]-constant
import math

class DoubleCB:
    def __call__(self, x, p):
        val = -99.
        t = float((x[0]-p[0])/p[1])
        if (t > -p[2] and t <  p[4]):
            val = math.exp(float(-0.5*t*t))
        elif (t <= -p[2] ):
            alpha1invn1 = p[2]/p[3]
            val = math.exp(float(-0.5*p[2]*p[2]))*math.pow(float(1.-alpha1invn1*(p[2]+t)),-p[3])
        elif (t >= p[4] ):
            alpha2invn2 = p[4]/p[5]
            val = math.exp(float(-0.5*p[4]*p[4]))*math.pow(float(1.-alpha2invn2*(p[4]-t)),-p[5])
        
        return p[6]*val;
