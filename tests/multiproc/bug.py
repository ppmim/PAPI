import multiprocessing, time           
def foo(x):                            
	time.sleep(3)                         

multiprocessing.Pool(1).apply(foo, [1])

