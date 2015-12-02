#ifndef SonicSculpturePiece_h
#define SonicSculpturePiece_h
#include <G3D/G3DAll.h>
#include "AudioSample.h"
class SonicSculpturePiece {
protected:
    CFrame          m_frame;
    Array<float>	m_audioSamples;
    Array<float>    m_deltas;
    Array<CFrame>   m_originalFrames;
    Array<CFrame>   m_transformedFrames;
    Array<float>    m_radii;

    bool            m_gpuUpdated;

    AttributeArray  m_vertices;
    AttributeArray  m_texCoord;
    AttributeArray  m_normals;
    IndexStream     m_indices;

    shared_ptr<UniversalMaterial> m_material;

	Sample sampleAudio(int windowIndex, float alpha);

    void getTransformedFramesFromOriginals();
    void uploadToGPU();
    void getCPUGeometry(CPUVertexArray& cpuVertexArray, Array<int>& cpuIndexArray) const;
    SonicSculpturePiece() : m_gpuUpdated(false) {}
public:
    shared_ptr<UniversalMaterial> material() const {
        return m_material;
    }

    void serialize(BinaryOutput& output) const;

    shared_ptr<AudioSample> getBaseAudioSample() const;

	shared_ptr<AudioSample> getAudioSampleFromRay(const Ray& ray);

    static shared_ptr<SonicSculpturePiece> create(shared_ptr<UniversalMaterial> material);

    void insert(const CFrame& frame, const float radius, const float delta, const Array<float>& samples = Array<float>());

    void render(RenderDevice* rd, const LightingEnvironment& env);

    void setShaderArgs(Args& args);

    static shared_ptr<SonicSculpturePiece> fromBinaryInput(shared_ptr<UniversalMaterial> material, BinaryInput& input);

    int size() const {
        return m_originalFrames.size();
    }
    
};


#endif