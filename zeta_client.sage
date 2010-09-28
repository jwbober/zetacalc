import os, sys
import pipes
import subprocess

number_of_threads = '2'

base_path = "/tmp/zeta_server"
new_client_pipe_name = os.path.join(base_path, "new_client_pipe")
return_pipe_name = os.path.join(base_path, "client_return_pipe")

ssh_base_command = "ssh -i /home/bober/key localhost "
ssh_new_client_command = ssh_base_command + " 'cat - > " + new_client_pipe_name + "'"
ssh_return_pipe_command = ssh_base_command + " 'cat " + return_pipe_name +"'"

new_client_pipeline = pipes.Template()
new_client_pipeline.append(ssh_new_client_command, '--')

return_pipeline = pipes.Template()
return_pipeline.append(ssh_return_pipe_command, '--')


def do_computation(work_unit):
    if work_unit=="shutdown":
        return False
    else:
        t, delta, N, start, length = work_unit
        print "doing computation from", start, "of length", length
        
        zeta_process = subprocess.Popen([ './zeta', str(t), str(delta), str(N), str(start), str(length), number_of_threads], stdout=subprocess.PIPE)

        line = zeta_process.stdout.readline()
        getting_result = False
        result = []
        while line != "":
            line = line.strip()
            if line == "BEGIN RESULT":
                getting_result = True
                line = zeta_process.stdout.readline()
                continue
            if getting_result:
                result.append(CC(line))
            else:
                print "subprocess output:", line
            line = zeta_process.stdout.readline()
        print result
        return vector(result)

def main():
    client_data = {'version': 1, 'return_pipe_name' : return_pipe_name, 'description' : 'test'}

    print "sending data to server"
    new_client_pipe = new_client_pipeline.open('/dev/null', 'w')
    new_client_pipe.write(dumps(client_data))
    new_client_pipe.close()

    print "waiting for response from server"
    #return_pipe = open(return_pipe_name, 'r')
    return_pipe = return_pipeline.open('/dev/null', 'r')
    client_id = loads(return_pipe.read())
    return_pipe.close()

    print "starting as client number", client_id
    input_pipe_name = os.path.join(base_path, "zeta_server_output" + str(client_id))
    output_pipe_name = os.path.join(base_path, "zeta_server_input" + str(client_id))
    
    input_pipeline = pipes.Template()
    input_pipeline.append(ssh_base_command + " 'cat " + input_pipe_name + "'", '--')

    output_pipeline = pipes.Template()
    output_pipeline.append(ssh_base_command + " 'cat - > " + output_pipe_name + "'", '--')

    while 1:
        output_pipe = output_pipeline.open('/dev/null', 'w')
        output_pipe.write("ready")
        output_pipe.close()
        #input_pipe = open(input_pipe_name, 'r')
        input_pipe = input_pipeline.open('/dev/null', 'r')
        message = input_pipe.read()
        input_pipe.close()

        result = do_computation(loads(message))
        if result is False:
            break
        else:
            output_pipe = output_pipeline.open('/dev/null', 'w')
            output_pipe.write(dumps(result))
            output_pipe.close()

            input_pipe = input_pipeline.open('/dev/null', 'r')
            message = input_pipe.read()
            input_pipe.close()
            if message != "thanks":
                print "problem"

            

        print "client number", client_id, "got message", message
        sleep(1)

if __name__== "__main__":
    main()
