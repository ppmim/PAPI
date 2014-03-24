"""
 PySQLite example for benchmarking
"""

ITERATIONS = 5

from pysqlite2 import dbapi2 as sqlite
import time

testdata = [(-7.345, 4200, "foooooooooooooooooooooxxxxxxxxxxxxxxx") for x in xrange(300000)]

con = None

def insert():
    global con
    con = sqlite.connect(":memory:")
    cur = con.cursor()
    cur.execute("create table foo(a,b,c)")
    cur.executemany("insert into foo(a,b,c) values (?,?,?)", testdata)

def fetch():
    cur = con.cursor()
    #cur.execute("select * from foo union select * from foo union select * from foo union select * from foo")
    cur.execute("select * from foo")
    #num = len(cur.fetchall())
    #for row in cur:
    #    print 'Row1=%f Row2=%d Row3=%s' % (row[0], row[1], row[2])
    rows=cur.fetchall()
    print "LEN=%d" %len(rows)
    print "A=%f" %rows[0].a

def timeit(func):
    t0 = time.time()
    func()
    return time.time() - t0

print "pysqlite", sqlite.version

sum = 0
for i in xrange(ITERATIONS):
    sum += timeit(insert)
print "average insert time: %f seconds" % (sum / ITERATIONS)

sum = 0
for i in xrange(ITERATIONS):
    sum += timeit(fetch)
print "average fetch time: %f seconds" % (sum / ITERATIONS)
print "-" * 50


