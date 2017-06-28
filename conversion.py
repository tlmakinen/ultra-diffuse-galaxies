# for converting Wittmann data to plot on Chandra evt's

def convertdistance(parsecdist):
    parsec_to_kpc = 0.332
    return (parsecdist / parsec_to_kpc) * 2

print(convertdistance(1.9))
