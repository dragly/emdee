#ifndef GLOGFALLBACK_H
#define GLOGFALLBACK_H

#include <iostream>

struct FallbackLogger {
    FallbackLogger() {
    }

    ~FallbackLogger() {
        std::cout << m_SS.str() << std::endl;
    }

public:
    // accepts just about anything
    template<class T>
    FallbackLogger &operator<<(const T &x) {
        m_SS << x;
        return *this;
    }
    // accepts endl
    typedef std::ostream& (manip)(std::ostream&);
    FallbackLogger &operator<<(manip& m) {
        (void)m;
        return *this;
    }
private:
    std::ostringstream m_SS;
};

#ifndef LOG
#define LOG(severity) FallbackLogger() << "LOG: "
#endif

#ifndef INFO
#define INFO
#endif

#ifndef DEBUG
#define DEBUG
#endif

#ifndef WARNING
#define WARNING
#endif

#ifndef ERROR
#define ERROR
#endif

#ifndef FATAL
#define FATAL
#endif

#endif // GLOGFALLBACK_H
