from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np 
import matplotlib.pyplot as plt

from . import logger

def calcSemis(totalMass, smallMass, largeMass, nSmall, 
    inner, outer, alpha, returnFigure=False):
    """
    calculate the position of bodies in a disk matching a given surface density profile

    the code assumes a bimodal mass distribution with small and large bodies
    commonly refered to as planetesimals and embryos

    Parameters
    ----------        
    totalMass : float
        the desired total mass in the disk in earth-masses
    smallMass, largeMass : float
        the mass of each small and large body in earth-masses
    nSmall : int
        the number of small bodies
    inner, outer : float
        the inner and outer orbital distance of the disk in AU
    alpha : float
        the disk surface density profile where
        sigma propto mass**-alpha
    returnFigure : bool, optional
        make a fugure showing the mass distribution

    Returns
    ----------   
    finalSemi: array
        the distribution of semimajor axes of the bodies
    finalMass: array
        the masses of the bodies at the given semimajor axes
    fig: figure object
        if the returnFigure keyword is True this is a matplotlib figure object
        depicting the 

    a warning will be printed if the final mass in bodies is more than 5 percent discrepent
    from the requested total mass
    """

    nLarge = int((totalMass - (nSmall * smallMass)) / largeMass)

    logger.info('{} large bodies'.format(nLarge))
    logger.info('{} total bodies'.format(nLarge+nSmall))

    missingMass = (totalMass - (nSmall * smallMass) - 
        (largeMass * nLarge))

    # produce a warning if there is more than 5% missing mass
    if missingMass >= totalMass * 0.05:
        logger.warning('{} missing mass'.format(missingMass))
    else:
        logger.info('{} missing mass'.format(missingMass))

    # how many times more massive is large relative to small
    massMultiplier = largeMass / smallMass

    # convert from surface density to radius
    power = alpha * 0.5

    span = outer - inner

    # calculate the initial positions using only small bodies
    unidist = np.linspace(0, 1, nSmall + ((largeMass * nLarge) / smallMass))
    powdist = unidist**(power)

    spacings = (powdist * span) + inner

    # we make the final disk by masking bodies
    maskSmall = np.ones_like(spacings, dtype=bool)
    maskLarge = np.zeros_like(spacings, dtype=bool)

    pos = len(unidist) / nLarge 
    for i in range(nLarge):
        maskLarge[int((i*pos)+(0.5*pos))] = True
        maskSmall[int((i*pos)+(0.5*pos)-(massMultiplier/2)):
            int((i*pos)+(0.5*pos)+(massMultiplier/2))] = False

    idx = np.argsort(np.r_[spacings[maskLarge],spacings[maskSmall]])
    finalSemi = np.r_[spacings[maskLarge],spacings[maskSmall]]
    finalMass = np.r_[np.ones_like(spacings[maskLarge]) * largeMass, 
        np.ones_like(spacings[maskSmall]) * smallMass]

    bodyType = np.where(finalMass > smallMass, 'EM', 'PL')

    if returnFigure:
        fig, [ax1,ax2]  = plt.subplots(2, 1, figsize=[9,7])
        ax1.scatter(np.arange(unidist.shape[0])[maskLarge],spacings[maskLarge],
            s=300,color='r',alpha=0.5,edgecolors='k')
        ax1.scatter(np.arange(unidist.shape[0])[maskSmall],spacings[maskSmall])
        ax1.grid()
        ax1.set_xlabel('Body number')
        ax1.set_ylabel('Semimajor axis')

        ax2.plot(finalSemi[idx],np.cumsum(finalMass[idx]), '.-')

        n = 1000
        npm = totalMass / n
        x = np.linspace(0,1,n)
        y = ((x ** (alpha * 0.5)) * (outer-inner)) + inner
        ax2.plot(y, np.cumsum(np.ones_like(y) * npm), color='r')
        ax2.set_xlabel('Cumulative mass')
        ax2.set_ylabel('Semimajor axis')


        return finalSemi[idx], finalMass[idx], bodyType[idx], fig


    return finalSemi[idx], finalMass[idx], bodyType[idx]


def calcMutualHill(a1, a2, m1, m2, mstar=1.0):
    """ calculate the mutual hill radius of to bodies

    Parameters
    ----------        
    a1, a2 : float
        semimajor axis of bodies m1 and m2
        in AU
    m1, m2 : float
        mass of bodies m1 and m2 in earth-masses
    mstar : float, optional
        mass of the central star in solar-masses

    Returns
    ----------   
    mutualHill: float
        the mutual hill radius
    """
    au2meter = 1.496E11
    mearth2kg = 5.9742E24
    a1_si = au2meter * a1
    a2_si = au2meter * a2 
    m1_si = mearth2kg * m1
    m2_si = mearth2kg * m2
    mstar_si = 1.9889E30 * mstar

    mutualHill = (0.5 * (a1_si + a2_si) * 
        ((m1_si + m2_si)/(3 * mstar_si))**(1/3)) / au2meter
    return mutualHill

def drawEIOoM(size=1,**kwargs):
    try:
        ecc = kwargs['ecc']
    except KeyError:
        ecc = np.random.uniform(0,0.01, size=size)

    try:
        inc = kwargs['inc']
    except KeyError:
        inc = np.random.uniform(0,0.5, size=size)

    try:
        littleOm = kwargs['littleOm']
    except KeyError:
        littleOm = np.random.uniform(0,360, size=size)

    try:
        bigOm = kwargs['bigOm']
    except KeyError:
        bigOm = np.random.uniform(0,360, size=size)

    try:
        meananom = kwargs['meananom']
    except KeyError:
        meananom = np.random.uniform(0,360, size=size)

    return ecc, inc, littleOm, bigOm, meananom




def writeHead():
    """write a header out in the format of a Mercury big.in file
    """
    
    headstr = """)O+_06 Big-body initial data  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------
style (Cartesian, Asteroidal, Cometary) = Asteroidal
epoch (in days) = 0.0
)---------------------------------------------------------------------
"""
    return headstr

def writeBody(bodyname, mass, semimajor, 
    ecc, inc, littleOm, bigOm, meananom):
    """write out a line for a single body in format of Mercury big.in file
    """
    bodystr= """ {} m={:.6e} r=0.1 d=3
  {:.6e} {:.7e} {:.4e} {:.4e} {:.4e} {:.4e} 0.0 0.0 0.0
""".format(bodyname, mass, semimajor, ecc, inc, littleOm, bigOm, meananom)
    return bodystr

def createBigin(totalMass, smallMass, largeMass, nSmall, 
    inner, outer, alpha, ):
    
    finalSemi, finalMass, bodyType = calcSemis(totalMass, smallMass, 
        largeMass, nSmall, inner, outer, alpha)

    nbodies = finalSemi.shape[0]
    if nbodies > 9999:
        logger.error('Too many bodies, need to change code in name bodies')
    ecc, inc, littleOm, bigOm, meananom = drawEIOoM(size=nbodies)

    outstr = writeHead()
    earthmass2sunmass = 3.003467e-6
    for i,a in enumerate(finalSemi):
        bodyname = bodyType[i] + '{:04}'.format(i)
        outstr += writeBody(bodyname, finalMass[i]*earthmass2sunmass, finalSemi[i], 
                        ecc[i], inc[i], littleOm[i], bigOm[i], meananom[i])
    #remove the final blank line
    outstr = outstr[:-1]
    return outstr





if __name__ == '__main__':

    ###
    totalmass = 5.0
    smallmass = 0.01
    largemass = 0.1
    nsmallbodies = 260
    inner = 0.2
    outer = 4.0
    alpha = 3/2
    #####

    finalSemi, finalMass, bodyType, fig = calcSemis(totalmass, smallmass, 
        largemass, nsmallbodies, inner, outer, alpha, returnFigure=True)
    fig.show()


