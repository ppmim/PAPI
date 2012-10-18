import multiprocessing, time           
def foo(t):                            
	time.sleep(t)                         
	return "Nada "

print "start"

r = multiprocessing.Pool(1).apply_async(foo, args=[3])
#r.wait()
print "Result=",r.get()

print "end !"
