#writing some useful math type functions
import operator

#like summation, but product (of a list of factors)
def prod(factors):
    return reduce(operator.mul, factors, 1)
