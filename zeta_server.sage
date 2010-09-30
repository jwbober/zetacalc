import os, sys
import threading

base_path = "/tmp/zeta_server"
new_client_pipe_name = os.path.join(base_path, "new_client_pipe")
return_pipe_name = os.path.join(base_path, "client_return_pipe")

class Computation:
    def __init__(self, t, delta, N):
        self.number_of_clients = 0
        self.clients = set([])
        self.t = t
        self.delta = delta
        self.N = N
        work_units = []
        self.end_point = floor(sqrt(t/(2 * pi)))
        self.work_unit_length = 10000000000
        self.next_index = 1
        self.lock = threading.Lock()

        self.result = vector([RR(0)] * N)

    def wait_for_new_clients(self):
        if not os.path.exists(new_client_pipe_name):
            os.mkfifo(new_client_pipe_name)
        if not os.path.exists(return_pipe_name):
            os.mkfifo(return_pipe_name)

        while 1:
            print "waiting for new client to connect"
            new_client_pipe = open(new_client_pipe_name, 'r')
            new_client_data = loads(new_client_pipe.read())
            new_client_pipe.close()
            self.create_new_client(new_client_data)

    def process_result(self, result):
        self.lock.acquire()
        self.result = self.result + result
        self.lock.release()

    def get_work_unit(self, ID):
        self.lock.acquire()
        if self.next_index > self.end_point:
            self.lock.release()
            self.unregister_client(ID)
            if len(self.clients) == 0:
                for x in self.result:
                    print x
                sys.exit()
            return "shutdown"
        start = self.next_index
        length = min(self.work_unit_length, self.end_point - start + 1)
        self.next_index = start + length
        self.lock.release()
        return (self.t, self.delta, self.N, start, length)

    def create_new_client(self, new_client_data):
        """
        new_client_data should be a dictionary, which should have
        an entry "version" and then other extra information which
        depends on the version.

        version 1:
            return_pipe_name (string): a filename. we will open this file
            and write to it to return some information to the client

            description (string): some descriptive string to describe this client
        """
        self.lock.acquire()

        print "creating new client"
        print new_client_data

        if(new_client_data['version'] == 1):
            new_client = Client(self.number_of_clients, new_client_data['description'], self)
            new_client.start()
            
            return_pipe = open(new_client_data['return_pipe_name'], 'w')
            return_pipe.write(dumps(self.number_of_clients))
            return_pipe.close()
            self.clients.add(number_of_clients)

            self.number_of_clients += 1
            self.clients.append(new_client)


        self.lock.release()
        
    def unregister_client(self, ID):
        self.lock.acquire()
        self.clients.remove(ID)
        self.lock.release()

class Client(threading.Thread):
    def __init__(self, ID, description, computation):
        threading.Thread.__init__(self)
        self.ID = ID
        self.desciption = description
        self.computation = computation

    def run(self):
        input_pipe_name = os.path.join(base_path, "zeta_server_input" + str(self.ID))
        output_pipe_name = os.path.join(base_path, "zeta_server_output" + str(self.ID))
        if not os.path.exists(input_pipe_name):
            os.mkfifo(input_pipe_name)
        if not os.path.exists(output_pipe_name):
            os.mkfifo(output_pipe_name)

        while 1:
            input_pipe = open(input_pipe_name, 'r')
            incoming = input_pipe.read()
            input_pipe.close()
            if incoming == "ready":
                print "giving work to client number", self.ID
                work_unit = self.computation.get_work_unit()
                output_pipe = open(output_pipe_name, 'w')
                output_pipe.write( dumps(work_unit) )
                output_pipe.close()
            elif incoming == "quitting":
                print "client number", self.ID, "shutting down."
                self.computation.unregister_client(self.ID)
            else:
                # assume that incoming is a vector containing the result of a computation
                result = loads(incoming)
                self.computation.process_result(result)
                output_pipe = open(output_pipe_name, 'w')
                output_pipe.write("thanks")
                output_pipe.close()

def main():
    if not os.path.exists(base_path):
        os.mkdir(base_path)
    
    c = Computation(100000000000000000000000, .01, 500)

    c.wait_for_new_clients()


if __name__ == "__main__":
     main()
