/**
 * Author: Yue Teng
 * Compile: g++ master.cpp -std=c++11 -o master
 * 
 * Description: The master program in charge of rgen, constructor, and 
 *              path_gener. It forks and execs all the programs, and 
 *              perform interprocess communication in the patern ---
 *              rgen -> constructor -> path_gener. STDIN dup2 
 *              pipe for path_gener.
 */

#include <iostream>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>

using namespace std;

int main(int argc, char *argv[])
{
    int fds_2a2[2];
    pipe(fds_2a2);
    int fds_rgen2a1[2];
    pipe(fds_rgen2a1);
    pid_t child_pid[4]; // fork 4 times

    // fork: parent and child (rgen)
    child_pid[0] = fork();
    if (child_pid[0] > 0)
    {
        // in parent
        // fork: parent and child (constructor)
        child_pid[1] = fork();
        if (child_pid[1] > 0)
        {
            // in parent
            // fork: parent and child (controller)
            child_pid[2] = fork();
            if (child_pid[2] > 0)
            {
                // in parent
                // fork: parent and child (path_gener)
                child_pid[3] = fork();
                if (child_pid[3] > 0)
                { // driver process
                    int wstatus;
                    pid_t cpid;
                    cpid = wait(&wstatus);
                    if (cpid > 0)
                    {
                        for (int i = 0; i < 4; ++i)
                        {
                            kill(child_pid[i], SIGKILL);
                        }
                        return 1;
                    }
                } // end driver
                else if (child_pid[3] == 0)
                { // processs ece650a2
                    close(fds_rgen2a1[0]);
                    close(fds_rgen2a1[1]);
                    close(fds_2a2[1]);
                    dup2(fds_2a2[0], STDIN_FILENO);
                    char const *arg_list[2] = {"./path_gener", NULL};
                    usleep(2); // sleep for 5 ms waiting for other processes
                    execvp(arg_list[0], (char *const *)arg_list);
                    
                    // finally
                    close(fds_2a2[0]);
                } // end ece650a2
                else
                { // error
                    return 1;
                } // end error
            }
            else if (child_pid[2] == 0)
            { // process inputctrl
                close(fds_rgen2a1[0]);
                close(fds_rgen2a1[1]);
                close(fds_2a2[0]);
                dup2(fds_2a2[1], STDOUT_FILENO);
                string msg; // input message e.g. s 2 10
                usleep(1);
                while (!cin.eof())
                {
                    getline(cin, msg);
                    if (cin.eof())
                    {
                        break;
                    }
                    fflush(stdin);
                    cout << msg << endl;
                }

                // finally
                close(fds_2a2[1]);
            } // end inputctrl
            else
            { // error
                return 1;
            } // end error
        }
        else if (child_pid[1] == 0)
        { // process ece650a1
            close(fds_2a2[0]);
            close(fds_rgen2a1[1]);
            dup2(fds_2a2[1], STDOUT_FILENO);
            dup2(fds_rgen2a1[0], STDIN_FILENO);
            char const *arg_list[3] = {"python3", "./constructor.py", NULL};
            usleep(1);
            execvp(arg_list[0], (char *const *)arg_list);

            // finally
            close(fds_rgen2a1[0]);
            close(fds_2a2[1]);
        } // end ece650a1
        else
        { // error
            return 1;
        } // end error
    }
    else if (child_pid[0] == 0)
    { // process rgen
        close(fds_2a2[0]);
        close(fds_2a2[1]);
        close(fds_rgen2a1[0]);
        dup2(fds_rgen2a1[1], STDOUT_FILENO);
        char const *arg_list[argc + 1];
        arg_list[0] = "./rgen";
        arg_list[argc] = NULL;
        for (int i = 1; i < argc; ++i)
        {
            arg_list[i] = argv[i];
        }
        execvp(arg_list[0], (char *const *)arg_list);

        // close the other end finally
        close(fds_rgen2a1[1]);
    } // end rgen
    else
    { // error
        return 1;
    } // end error
    return 0;
}
