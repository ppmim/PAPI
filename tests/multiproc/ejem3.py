from multiprocessing import Process, Pipe, Pool


def f(conn):
    conn.send([42, None, 'hello'])
    conn.close()

if __name__ == '__main__':
    parent_conn, child_conn = Pipe()
    p = Process(target=f, args=(child_conn,))
    p.start()
    p.join()
    print parent_conn.recv()   # prints "[42, None, 'hello']"
   
