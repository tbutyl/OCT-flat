import multiprocessing as mp
#from os import getpid


def cube(x):
    #(getpid())
    return x**3

def runner(i):
    
    with mp.Pool() as pool:
        result = pool.map(cube, range(i))
    
    return result

if __name__ == '__main__':

    result = runner(100)
    print(result)


#results = [cube(i) for i in range(7)]
#print(results)