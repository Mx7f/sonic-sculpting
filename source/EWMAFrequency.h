#ifndef EWMAFrequency_h
#define EWMAFrequency_h

#include <G3D/G3DAll.h>
/** Exponentially-weighted moving average of frequencies */
struct EWMAFrequency {
    Array<float>        cpuData;
    shared_ptr<Texture> gpuData;
    // Update rule: freq_ewma = lerp(freq_current, freq_ewma, alpha)
    float               alpha;

    void init(float a, const Array<float>& newData, const String& name) {
        alpha = a;
        cpuData.appendPOD(newData);
        gpuData = Texture::createEmpty(name, newData.size(), 1, ImageFormat::R32F());
        upload();
    }

    void update(const Array<float>& newData) {
        alwaysAssertM(newData.size() == cpuData.size(), "Must have same size data for EWMAFrequency update");
        for (int i = 0; i < newData.size(); ++i) {
            cpuData[i] = lerp(newData[i], cpuData[i], alpha);
        }
        upload();
    }

    void upload() {
        debugAssert(notNull(gpuData));
        shared_ptr<CPUPixelTransferBuffer> ptb = CPUPixelTransferBuffer::fromData(cpuData.size(), 1, ImageFormat::R32F(), cpuData.getCArray());
        gpuData->update(ptb);
    }
};
#endif
