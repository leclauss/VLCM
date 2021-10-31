#ifndef VLCM_EXCEPTION_H
#define VLCM_EXCEPTION_H

#include <exception>
#include <string>

class VLCMException : public std::exception {
public:
    explicit VLCMException(std::string message) : _msg(std::move(message)) {}

    ~VLCMException() override = default;

    VLCMException(const VLCMException &copyFrom) = default;

    VLCMException &operator=(const VLCMException &copyFrom) = default;

    VLCMException(VLCMException &&) = default;

    VLCMException &operator=(VLCMException &&) = default;

    const char *what() const noexcept override { return _msg.c_str(); }

protected:
    std::string _msg;
};


#endif //VLCM_EXCEPTION_H
