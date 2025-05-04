#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <atomic>
#include <chrono>

class SpinProgress {
public:
    enum SpinType {
        Classic,
        Braille,
        Dots,
        Circle,
        Arrow,
        Cross,
        Blocks,
        CheckMark,
        Unicode
    };

    SpinProgress(SpinType type = Classic)
        : is_completed(false), custom_message(""), spinner_thread(&SpinProgress::spin, this) {
        set_spin_chars(type);
    }

    ~SpinProgress() {
        if (spinner_thread.joinable()) {
            complete();
            spinner_thread.join();
        }
    }

    void set_spin_chars(SpinType type) {
        switch (type) {
        case Classic:
            spin_chars = { "|", "/", "-", "\\" };
            break;
        case Dots:
            spin_chars = { ".", "o", "O", "*", "@" };
            break;
        
        }
    }

    void set_message(const std::string& message) {
        custom_message = message;
    }

    void complete(const std::string& message = "Completed") {
        is_completed.store(true);
        if (spinner_thread.joinable()) {
            spinner_thread.join();
        }
        std::cout << "\r\033[K";  // カーソルを行の先頭に戻し、行をクリア
        std::cout << message << std::endl;
    }

private:
    std::vector<std::string> spin_chars;  // std::vector<char> から std::vector<std::string> に変更
    std::atomic<bool> is_completed;
    std::string custom_message;
    std::thread spinner_thread;

    void spin() {
        size_t spin_index = 0;
        while (!is_completed.load()) {
            std::cout << "\r" << spin_chars[spin_index % spin_chars.size()] << " " << custom_message << std::flush;
            spin_index++;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
};