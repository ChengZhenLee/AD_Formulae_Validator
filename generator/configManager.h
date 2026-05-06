#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

#include <string>
#include <map>
#include <vector>


class ConfigManager{
    private:
        std::map<std::string, std::string> config;
        ConfigManager() {};

    public:
        static ConfigManager& getInstance() {
            static ConfigManager instance;
            return instance;
        }

        void load(const std::string &filename);
        std::string getType();
        size_t getTangentShape();
        size_t getAdjointShape();
        size_t getXShape();
        size_t getYShape();
        std::string getPrimalFunctionName();
};

#endif