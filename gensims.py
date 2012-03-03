import random

nsims = 50
fname = 'run{}.sh'
cmd = '''
#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
$HOME/ramseymc/demon.out %d %d %d 1000 %d > 4-6-35_%d-%d-%d_%d.log\n
'''
base = [1000, 100000, 20]

for i in xrange(nsims):
    f = open(fname.format(i), 'w')
    seed = random.randint(0, 2**31-1)
    params = []
    for param in base:
        r = random.random()
        if r < 0.167:
            params.append(param*1.5)
        elif r < 0.333:
            params.append(param*0.667)
        else:
            params.append(param)

    f.write(cmd % (2*tuple(params+[seed])))
    f.close()





