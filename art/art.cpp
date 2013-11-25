#include "art.h"
using namespace std;

extern "C" {
    art* create(reconVars* pVars,
                float *pAngle,
                float *pCenter,
                float *pInput,
                float *pOutput) {
        reconVars pVars_;
        float *pAngle_ = NULL;
        float *pCenter_ = NULL;
        float *pInput_ = NULL;
        //float *pOutput_ = NULL;

        int bytesAngle = (pVars -> numProjections)*sizeof(float);
        int bytesCenter = (pVars -> numProjections)*sizeof(float);
        int bytesInput = (pVars -> numProjections)*(pVars -> numSlices)*(pVars -> numPixels)*sizeof(float);
        //int bytesOutput = (pVars -> numSlices)*(pVars -> imgN)*(pVars -> imgN)*sizeof(float);

        pAngle_ = (float *)malloc(bytesAngle);
        pCenter_ = (float *)malloc(bytesCenter);
        pInput_ = (float *)malloc(bytesInput);
        //pOutput_ = (float *)malloc(bytesOutput);

        memcpy(&pVars_, pVars, sizeof(pVars_));
        memcpy(pAngle_, pAngle, bytesAngle);
        memcpy(pCenter_, pCenter, bytesCenter);
        memcpy(pInput_, pInput, bytesInput);
        //memcpy(pOutput, pOutput, bytesOutput);

        return new art(&pVars_,
                       pAngle_,
                       pCenter_,
                       pInput_,
                       pOutput);
        }

    void reconstruct(art *art) {
        art -> reconstruct();
        }
    } // extern "C"

art::art(reconVars* pVars,
         float *pAngle,
         float *pCenter,
         float *pInput,
         float *pOutput)
  : vpVars(pVars),
    vNumPixels(pVars -> numPixels),
    vNumSlices(pVars -> numSlices),
    vNumProjections(pVars -> numProjections),
    vImgN(pVars -> imgN),
    vImgLen(pVars -> imgLen),
    vDetLen(pVars -> detLen),
    vpAngle(pAngle),
    vpCenter(pCenter),
    vpInput(pInput),
    vpOutput(pOutput) {
    }

void art::reconstruct() {
    int m, n, k, l;
    float imgGridSize;
    float xi[vImgN+1];
    float yi[vImgN+1];
    float detx0[vNumPixels];
    float dety0[vNumPixels];
    float srcx0[vNumPixels];
    float srcy0[vNumPixels];
    float detx;
    float dety;
    float srcx;
    float srcy;
    std::vector<float> ax, ay, alpha;
    std::vector<float> axCleaned, ayCleaned;
    float alphaDiff;
    float amin;
    float amax;
    std::vector<float> xk, yk;
    std::vector<float> x0, y0;
    std::vector<float> dist;
    float dist2;
    int row, col;
    int indIn;
    std::vector<int> indOut;
    float distvpOutput;

    // Convert angles in radians.
    for (m = 0; m < vNumProjections; m++) {
        vpAngle[m] *= PI/180;
        }

    //for (m = 0; m < vNumPixels*vNumSlices*vNumProjections; m++) {
    //    std::cout << m << std::endl;
    //    std::cout << vpInput[m] << std::endl;
    //    }

    // Image grid size
    imgGridSize = float(vNumPixels) / vImgN;

    // Find grid intersections.
    for (m = 0; m < vImgN+1; m++) {
        xi[m] = -float(vNumPixels)/2 + m*imgGridSize;
        yi[m] = -float(vNumPixels)/2 + m*imgGridSize;
        }

    // Initial source/detector pair locations.
    for (m = 0; m < vNumPixels; m++) {
        detx0[m] = -1e4;
        dety0[m] = -float(vNumPixels-1)/2 + m;
        srcx0[m] = 1e4;
        srcy0[m] = -float(vNumPixels-1)/2 + m;
        }
    //std::cout << "srcx0: " << std::endl;
    //for (k = 0; k < vNumPixels; k++) {
    //    std::cout << srcx0[k] << std::endl;
    //    }
    //std::cout << "srcy0: " << std::endl;
    //for (k = 0; k < vNumPixels; k++) {
    //    std::cout << srcy0[k] << std::endl;
    //    }
    //std::cout << "detx0: " << std::endl;
    //for (k = 0; k < vNumPixels; k++) {
    //    std::cout << detx0[k] << std::endl;
    //    }
    //std::cout << "dety0: " << std::endl;
    //for (k = 0; k < vNumPixels; k++) {
    //    std::cout << dety0[k] << std::endl;
    //    }
    //std::cout << vNumProjections << std::endl;
    //for (k = 0; k < vNumProjections; k++) {
    //    std::cout << vpAngle[k] << std::endl;
    //    }
    //std::cout << "...." << std::endl;

    for (m = 0; m < vNumProjections; m++) { //vNumProjections
        std::cout << m << std::endl;
        // Update detector coordinates.
        for (n = 0; n < vNumPixels; n++) { //vNumPixels
            //std::cout << n << std::endl;
            indIn = m + (n * vNumProjections);
            //std::cout << indIn << std::endl;
            // Init vectors.
            if (!ax.empty()) {
                ax.clear();
                }
            if (!ay.empty()) {
                ay.clear();
                }
            if (!axCleaned.empty()) {
                axCleaned.clear();
                }
            if (!ayCleaned.empty()) {
                ayCleaned.clear();
                }
            if (!alpha.empty()) {
                alpha.clear();
                }
            if (!xk.empty()) {
                xk.clear();
                }
            if (!yk.empty()) {
                yk.clear();
                }
            if (!dist.empty()) {
                dist.clear();
                }
            if (!x0.empty()) {
                x0.clear();
                }
            if (!y0.empty()) {
                y0.clear();
                }
            if (!indOut.empty()) {
                indOut.clear();
                }

            //std::cout << n << std::endl;
            srcx = srcx0[n] * cos(vpAngle[m]) + srcy0[n] * sin(vpAngle[m]);
            srcy = srcx0[n] * sin(vpAngle[m]) - srcy0[n] * cos(vpAngle[m]);
            detx = detx0[n] * cos(vpAngle[m]) + dety0[n] * sin(vpAngle[m]);
            dety = detx0[n] * sin(vpAngle[m]) - dety0[n] * cos(vpAngle[m]);
            //std::cout << srcx << std::endl;
            //std::cout << srcy << std::endl;
            //std::cout << detx << std::endl;
            //std::cout << dety << std::endl;

            // Calculate the alpha values.
            for (k = 0; k < vImgN+1; k++) {
                //std::cout << k << std::endl;
                ax.push_back((xi[k] - srcx) / (detx - srcx));
                ay.push_back((yi[k] - srcy) / (dety - srcy));
                }
            //std::cout << "ax size: " << std::endl;
            //std::cout << ax.size() << std::endl;
            //std::cout << "ax: " << std::endl;
            //for (k = 0; k < ax.size(); k++) {
            //    std::cout << ax[k] << std::endl;
            //    }
            //std::cout << "ay size: " << std::endl;
            //std::cout << ay.size() << std::endl;
            //std::cout << "ay: " << std::endl;
            //for (k = 0; k < ay.size(); k++) {
            //    std::cout << ay[k] << std::endl;
            //    }

            amin = fmaxf(0, fmaxf(fminf(ax[0], ax[vImgN]), fminf(ay[0], ay[vImgN])));
            amax = fminf(1, fminf(fmaxf(ax[0], ax[vImgN]), fmaxf(ay[0], ay[vImgN])));

            for (k = 0; k < vImgN+1; k++) {
                if ((ax[k] > amin) && (ax[k] < amax)) {
                    axCleaned.push_back(ax[k]);
                    }
                if ((ay[k] > amin) && (ay[k] < amax)) {
                    ayCleaned.push_back(ay[k]);
                    }
                }
            //std::cout << "axCleaned size: " << std::endl;
            //std::cout << axCleaned.size() << std::endl;
            //std::cout << "axCleaned: " << std::endl;
            //for (k = 0; k < axCleaned.size(); k++) {
            //    std::cout << axCleaned[k] << std::endl;
            //    }
            //std::cout << "ayCleaned size: " << std::endl;
            //std::cout << ayCleaned.size() << std::endl;
            //std::cout << "ayCleaned: " << std::endl;
            //for (k = 0; k < ayCleaned.size(); k++) {
            //    std::cout << ayCleaned[k] << std::endl;
            //    }

            std::merge(axCleaned.begin(), axCleaned.end(),
                       ayCleaned.begin(), ayCleaned.end(),
                       std::back_inserter(alpha));
            std::sort(alpha.begin(), alpha.end());
            //std::cout << "alpha size: " << std::endl;
            //std::cout << alpha.size() << std::endl;
            //std::cout << "alpha: " << std::endl;
            //for (k = 0; k < alpha.size(); k++) {
            //    std::cout << alpha[k] << std::endl;
            //    }

            // Inside ROI do:
            if (alpha.size() > 1) {
                // Calculate intersections.
                for (k = 0; k < alpha.size(); k++) {
                    xk.push_back(srcx + alpha[k] * (detx - srcx));
                    yk.push_back(srcy + alpha[k] * (dety - srcy));
                    }
                //std::cout << "xk.size(): " << std::endl;
                //std::cout << xk.size() << std::endl;
                //std::cout << "xk: " << std::endl;
                //for (k = 0; k < alpha.size(); k++) {
                //    std::cout << xk[k] << std::endl;
                //    }
                //std::cout << "yk: " << std::endl;
                //for (k = 0; k < alpha.size(); k++) {
                //    std::cout << yk[k] << std::endl;
                //    }

                // Calculate lengths.
                for (k = 0; k < alpha.size()-1; k++) {
                    dist.push_back(sqrt(pow(xk[k+1] - xk[k], 2) + pow(yk[k+1] - yk[k], 2)));
                    }
                //std::cout << "dist: " << std::endl;
                //for (k = 0; k < dist.size(); k++) {
                //    std::cout << dist[k] << std::endl;
                //    }

                // Calculate length^2.
                dist2 = 0;
                for (k = 0; k < dist.size(); k++) {
                    dist2 += dist[k] * dist[k];
                    }
                //std::cout << "dist2: " << std::endl;
                //std::cout << dist2 << std::endl;

                for (k = 0; k < dist.size(); k++) {
                    // Take the mid point of the crossings to
                    // find the line segment lying on that pixel.
                    x0.push_back((xk[k] + xk[k+1]) / 2);
                    y0.push_back((yk[k] + yk[k+1]) / 2);

                    // Update image based on ART.
                    row = floor(x0[k] / imgGridSize + float(vImgN)/2);
                    col = floor(y0[k] / imgGridSize + float(vImgN)/2);
                    indOut.push_back(col + (row * vImgN));
                    //std::cout << float(vImgN)/2 << std::endl;
                    }
                //std::cout << "indOut: " << std::endl;
                //for (k = 0; k < indOut.size(); k++) {
                //    std::cout << indOut[k] << std::endl;
                //    }

                distvpOutput = 0;
                for (k = 0; k < dist.size(); k++) {
                    //std::cout << indOut[k] << std::endl;
                    distvpOutput += vpOutput[indOut[k]] * dist[k];
                    }
                //std::cout << "distvpOutput:" << std::endl;
                //std::cout << distvpOutput << std::endl;

                for (k = 0; k < dist.size(); k++) {
                    vpOutput[indOut[k]] += (vpInput[indIn] - distvpOutput) / dist2 *  dist[k];
                    }

                    //std::cout << "::" << std::endl;
                    //std::cout << indOut << std::endl;
                    //std::cout << ((vpInput[indIn] - dist[k] * vpOutput[indOut]) / dist2 *  dist[k]) << std::endl;

                    //for (l = 0; l < vNumSlices; l++) {
                    //    indIn = m + (n - 1) * vImgN + l * (vNumPixels * vNumPixels);
                    //    indOut = (row + (col - 1) * vImgN) + l * (vImgN * vImgN);
                    //    vpOutput[indOut] += (vpInput[indIn] - dist * vpOutput[indOut]) / (dist * dist) * dist;
                    //    }
                }

            }
        }
    }
