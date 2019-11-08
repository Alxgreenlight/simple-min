#ifndef HELPER_ALX_HPP
#define HELPER_ALX_HPP
#include <iostream>

#define benchlib 0
#define ltool 1

namespace helper{
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