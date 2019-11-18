#ifndef HELPER_ALX_HPP
#define HELPER_ALX_HPP
#include <iostream>

#define benchlib 0
#define ltool 1

namespace helper{

    static std::string bhelp = "This program should be used with parameter:\nHere you can set accuracy with which calculation\
should be runned. Just type foolowing:\n./program_name eps\nWhere eps is float, e.g. 0.01\n\
Now let's restart program with right parameter";
    static std::string lthelp = "This program should be used with parameter:\n\
It depends on the benchmark type you want to run. Just type:\n\
./program_name [gN|ml]\n\
Where gN means GKLS with dimensionality N (integer, e.g. g2), and ml means mathexplib (used functions defined in code)\n\
Now let's restart program with right parameter";

    void help(int tool){
        switch (tool){
            case benchlib:
                std::cout << bhelp << std::endl;
                break;
            case ltool:
                std::cout << lthelp << std::endl;
                break;
            default:
                std::cout << std::endl;
        }
    }

    void progress_bar(int done, int todo){
        int barWidth = 70;
        float progress = 1.0*done/todo;
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();

        if (done == todo){
            std::cout << std::endl;
        }
        return;
    }
}

#endif