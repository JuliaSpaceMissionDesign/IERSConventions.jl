# TODO: improve documentation of this types! 

export  iers1996, 
        iers2003a, iers2003b, 
        iers2010a, iers2010b, CPNc, CPNd

abstract type IERSModel end 

# IERS 1996
# ========================

struct IERS1996 <: IERSModel end 

""" 
    iers1996

The singleton instance of type `IERS1996`, representing the IERS 1996 family of models.
"""
const iers1996 = IERS1996() 


# IERS 2003 
# ========================

abstract type IERS2003 <: IERSModel end 

struct IERS2003A <: IERS2003 end

struct IERS2003B <: IERS2003 end 

"""
    iers2003a

The singleton instance of type `IERS2003A`, representing the IERS 2003A family of models.
"""
const iers2003a = IERS2003A()

"""
    iers2003b

The singleton instance of type `IERS2003B`, representing the IERS 2003B family of models.
"""
const iers2003b = IERS2003B()


# IERS 2010
# ========================

abstract type IERS2010 <: IERSModel end 

struct IERS2010A <: IERS2010 end 

struct IERS2010B <: IERS2010 end 

abstract type CPNModel <: IERS2010 end 

struct CPNC <: CPNModel end 

struct CPND <: CPNModel end 


"""
    iers2010a

The singleton instance of type `IERS2010A`, representing the IERS 2010A family of models.
"""
const iers2010a = IERS2010A()

"""
    iers2010b

The singleton instance of type `IERS2010B`, representing the IERS 2010B family of models.

!!! note 
    This is not an official IERS model.
"""
const iers2010b = IERS2010B()

"""
    CPNc

The singleton instance of type `CPNC`, representing the concise CPNc from Capitaine & 
Wallace, Concise CIO based precession-nutation formulations, (2008). This model truncates 
the X, Y series to deliver an accuracy of few mas.

!!! note 
    This is not an official IERS model.
"""
const CPNc = CPNC()

"""
    CPNd

The singleton instance of type `CPND`, representing the concise CPNd from Capitaine & 
Wallace, Concise CIO based precession-nutation formulations, (2008). This model truncates 
the X, Y series to deliver an accuracy of few arcseconds.

!!! note 
    This is not an official IERS model.
"""
const CPNd = CPND()


# This is a type alias to group the highest accuracy models
IERSAModels = Union{IERS2003A, IERS2010A}

# This is a type alias to group the approximate models originating from IAU 2000B 
IERSBModels = Union{IERS2003B, IERS2010B}