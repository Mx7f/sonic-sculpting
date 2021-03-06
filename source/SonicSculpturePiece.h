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

    Vector3         m_generalDirection;
    bool            m_frozen;

    bool            m_gpuUpdated;

    AttributeArray  m_vertices;
    AttributeArray  m_texCoord;
    AttributeArray  m_normals;
    AttributeArray  m_boneIndices;
    AttributeArray  m_boneWeights;
    IndexStream     m_indices;

    shared_ptr<Texture> m_boneTexture;

    shared_ptr<UniversalMaterial> m_material;

    void uploadToGPU();
    void getCPUGeometry(CPUVertexArray& cpuVertexArray, Array<int>& cpuIndexArray) const;
    SonicSculpturePiece() : m_gpuUpdated(false) {}
public:
    shared_ptr<UniversalMaterial> material() const {
        return m_material;
    }

    void onSimulation(RealTime rdt, SimTime sdt, SimTime idt);

    Any toAny(const String& binaryFilename) const;

    static shared_ptr<SonicSculpturePiece> create(shared_ptr<UniversalMaterial> material, const Any& any);

    shared_ptr<AudioSample> getBaseAudioSample() const;

	shared_ptr<AudioSample> getAudioSampleFromRay(const Ray& ray);

    float minValueAlongRay(const Ray& ray) const;

    void minMaxValue(const CFrame& frame, Point3& minP, Point3& maxP) const;

    static shared_ptr<SonicSculpturePiece> create(shared_ptr<UniversalMaterial> material);

    void insert(const CFrame& frame, const float radius, const float delta, const Array<float>& samples = Array<float>());

    void render(RenderDevice* rd, const LightingEnvironment& env);

    void setShaderArgs(Args& args);



    int size() const {
        return m_originalFrames.size();
    }

    void setFrozen(bool b) {
        m_frozen = b;
    }
    
};


#endif