#ifndef PROGRESSREPORT_H
#define PROGRESSREPORT_H

#include <exception>

#include <stdio.h>
#include <curl/curl.h>
#include <time.h>
#include <iostream>
#include <string>
#include <thread>

using namespace std;

class ProgressReporter {
public:
    ProgressReporter(string username, string runName);

    void reportProgress(double progress);
    static void callAndCleanCurl(CURL* curl);
private:
    double m_lastReportTime;
    string m_username;
    int m_runID;
    string m_runName;
};

inline ProgressReporter::ProgressReporter(string username, string runName) :
    m_lastReportTime(0),
    m_username(username),
    m_runName(runName)
{
    srand(time(NULL));
    m_runID = rand();
}

inline void ProgressReporter::reportProgress(double progress)
{
    time_t timer;
    struct tm someTime;
    double seconds;

    someTime.tm_hour = 0;   someTime.tm_min = 0; someTime.tm_sec = 0;
    someTime.tm_year = 100; someTime.tm_mon = 0; someTime.tm_mday = 1;
    time(&timer);  /* get current time; same as: timer = time(NULL)  */

    seconds = difftime(timer,mktime(&someTime));
    if(progress != 0 && progress != 1 && seconds - 10 < m_lastReportTime) { // skip all reports that come too soon, except the first and the last
        return;
    }
    if(progress < 0 || progress > 1) {
        cerr << "Progress must be between 0 and 1." << endl;
        throw new std::exception();
    }
    CURL *curl;
    curl = curl_easy_init();
    if(curl) {
        stringstream url;
        url << "http://compphys.dragly.org/wp-content/plugins/run-reporter/submit.php?";
        url << "user=" << m_username;
        url << "&runName=" << m_runName;
        url << "&runID=" << m_runID;
        url << "&progress=" << progress;
        curl_easy_setopt(curl, CURLOPT_URL, url.str().c_str());
        curl_easy_setopt(curl, CURLOPT_TIMEOUT_MS, 5000);
        /* example.com is redirected, so we tell libcurl to follow redirection */
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);

        /* Perform the request, res will get the return code */
        //        res = curl_easy_perform(curl);
        /* Check for errors */
        //        if(res != CURLE_OK)
        //            fprintf(stderr, "curl_easy_perform() failed: %s\n",
        //                    curl_easy_strerror(res));

        thread test(ProgressReporter::callAndCleanCurl, curl);
        if(progress == 1) {
            test.join();
        } else {
            test.detach();
        }
    }
    m_lastReportTime = seconds;
}

inline void ProgressReporter::callAndCleanCurl(CURL *curl)
{
    /* Perform the request, res will get the return code */
    CURLcode res;
    res = curl_easy_perform(curl);
    /* Check for errors */
    if(res != CURLE_OK) {
        fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
    }
    /* always cleanup */
    curl_easy_cleanup(curl);
}
#endif // PROGRESSREPORT_H
