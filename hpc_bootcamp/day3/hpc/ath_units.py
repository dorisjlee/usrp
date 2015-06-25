
class Units_LMT(object):
    def __init__(self, Lcode=1.0,Mcode=1.0,Tcode=1.0):
        super(Units_LMT, self).__init__()
        cm        = 1.0/Lcode
        g         = 1.0/Mcode
        s         = 1.0/Tcode
        dyne      = g*cm/s**2
        erg       = g*cm**2/s**2
        K         = 1.0
        self.cm   = cm
        self.g    = g
        self.s    = s
        self.dyne = dyne
        self.erg  = erg
        self.K    = K
  
        #  Convert common physical constants/units from cgs to new units
        self.G    = 6.67259e-8 * cm**3/(g*s**2);
        self.Msun = 1.9891e+33 * g;
        self.Lsun = 3.8268e+33 * erg/s;
        self.Myr  = 3.155815e+13 * s;
        self.pc   = 3.085678e+18 * cm;
        self.kpc  = 3.085678e+21 * cm;
        self.kms  = 1.0e+5 * cm/s;
        self.mH   = 1.6733e-24 * g;
        self.aR   = 7.5646e-15 * erg/(cm**3*K**4);
        self.kB   = 1.380658e-16 * erg/K;
        self.c    = 2.99792458e+10 * cm/s;
        self.NA   = 6.0221367e23;
        
class Units_LDVK(object):
    def __init__(self, Lcode=1.0,Dcode=1.0,Vcode=1.0,Kcode=1.0):
        super(Units_LDVK, self).__init__()
        cm        = 1.0/Lcode
        g         = 1.0/(Dcode*Lcode**3)
        s         = Vcode/Lcode
        dyne      = g*cm/s**2
        erg       = g*cm**2/s**2
        K         = 1.0/Kcode
        self.cm   = cm
        self.g    = g
        self.s    = s
        self.dyne = dyne
        self.erg  = erg
        self.K    = K
  
        #  Convert common physical constants/units from cgs to new units
        self.G    = 6.67259e-8 * cm**3/(g*s**2);
        self.Msun = 1.9891e+33 * g;
        self.Lsun = 3.8268e+33 * erg/s;
        self.Myr  = 3.155815e+13 * s;
        self.pc   = 3.085678e+18 * cm;
        self.kpc  = 3.085678e+21 * cm;
        self.kms  = 1.0e+5 * cm/s;
        self.mH   = 1.6733e-24 * g;
        self.aR   = 7.5646e-15 * erg/(cm**3*K**4);
        self.kB   = 1.380658e-16 * erg/K;
        self.c    = 2.99792458e+10 * cm/s;
        self.NA   = 6.0221367e23;