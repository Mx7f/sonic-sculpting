#include "SonicSculpturePiece.h"

void SonicSculpturePiece::getTransformedFramesFromOriginals() {
    m_transformedFrames.fastClear();
    m_transformedFrames.resize(m_originalFrames.size());
    float totalDelta = 0.0f;
    // TODO: Find more elegant algorithm
    for (int i = m_originalFrames.size() - 1; i >= 0; --i) {
        m_transformedFrames[i] = m_originalFrames[i];
        m_transformedFrames[i].translation += m_originalFrames[i].lookVector() * totalDelta;
        totalDelta += m_deltas[i];
    }
}

void SonicSculpturePiece::getCPUGeometry(CPUVertexArray & cpuVertexArray, Array<int>& cpuIndexArray) const {
    int sideCount = 6;
    int segmentCount = m_transformedFrames.size();
    int currentPointIndex = 0;
    cpuVertexArray.vertex.resize(segmentCount*sideCount);
    cpuVertexArray.hasBones = true;
    cpuVertexArray.boneIndices.resize(segmentCount*sideCount);
    cpuVertexArray.boneWeights.resize(segmentCount*sideCount);

    for (int j = 0; j < segmentCount; ++j) {
        const CFrame& frame = m_transformedFrames[j];
        const float thickness = m_radii[j];
        const Point3& currentPoint = Point3(0.0, 0.0, 0.0);

        for (int k = 0; k < sideCount; ++k) {
            int index = currentPointIndex * sideCount + k;
            CPUVertexArray::Vertex& v = cpuVertexArray.vertex[index];
            float theta = k * 2.0f * pif() / sideCount;
            v.normal = Vector3(cos(theta), sin(theta), 0.0f);
            v.position = currentPoint + v.normal * thickness;
            v.texCoord0 = Point2(float(j) / float(segmentCount-1), float(k) / (sideCount-1));
            v.tangent = Vector4(Vector3::unitZ(), 1.0f);

            cpuVertexArray.boneIndices[index] = Vector4int32(j, j, j, j);
            cpuVertexArray.boneWeights[index] = Vector4(1.0f, 0.0f, 0.0f, 0.0f);

        }

        if (j > 0) {
            int currStart = currentPointIndex * sideCount;
            int prevStart = currStart - sideCount;
            for (int k = 0; k < sideCount; ++k) {
                const int next = (k + 1) % sideCount;
                cpuIndexArray.append(prevStart + k, prevStart + next, currStart + k);
                cpuIndexArray.append(prevStart + next, currStart + next, currStart + k);
            }
        }

        ++currentPointIndex;
    }
}



static void uploadBones
(const shared_ptr<Texture>& boneTexture,
    const Array<CFrame>&    boneMatrices) {

    if (notNull(boneTexture)) {
        // Copy Bones to GPU
        const shared_ptr<CPUPixelTransferBuffer>& pixelBuffer =
            CPUPixelTransferBuffer::create(boneTexture->width(),
                boneTexture->height(),
                boneTexture->format());
        Vector4* row0 = (Vector4*)pixelBuffer->row(0);
        Vector4* row1 = (Vector4*)pixelBuffer->row(1);

        for (int i = 0; i < boneMatrices.size(); ++i) {
            const CFrame& boneFrame = boneMatrices[i];
            /* Unoptimized but readable version:
            const Matrix4& boneMatrix = boneFrame.toMatrix4();
            *row0   = boneMatrix.column(0);
            *row1   = boneMatrix.column(1);
            ++row0; ++row1;
            *row0   = boneMatrix.column(2);
            *row1   = boneMatrix.column(3);
            ++row0; ++row1;
            */
            const Matrix3& R = boneFrame.rotation;
            const Vector3& T = boneFrame.translation;
            row0->x = R[0][0];
            row0->y = R[1][0];
            row0->z = R[2][0];
            row0->w = 0.0f;
            row1->x = R[0][1];
            row1->y = R[1][1];
            row1->z = R[2][1];
            row1->w = 0.0f;
            ++row0; ++row1;

            row0->x = R[0][2];
            row0->y = R[1][2];
            row0->z = R[2][2];
            row0->w = 0.0f;
            row1->x = T.x;
            row1->y = T.y;
            row1->z = T.z;
            row1->w = 1.0f;
            ++row0; ++row1;
        }

        boneTexture->update(pixelBuffer);
    }
}

void SonicSculpturePiece::uploadToGPU() {
    getTransformedFramesFromOriginals();

    CPUVertexArray cpuVertexArray;
    Array<int> cpuIndexArray;
    getCPUGeometry(cpuVertexArray, cpuIndexArray);
    AttributeArray ignore;
    cpuVertexArray.copyToGPU(m_vertices, m_normals, ignore, m_texCoord, ignore, ignore, m_boneIndices, m_boneWeights, VertexBuffer::WRITE_ONCE);
    m_indices = IndexStream(cpuIndexArray, VertexBuffer::create(sizeof(int)*cpuIndexArray.size() + 8));
    m_boneTexture = Texture::createEmpty("Sonic Sculpture Bone Texture", m_transformedFrames.size() * 2, 2, ImageFormat::RGBA32F());

    uploadBones(m_boneTexture, m_transformedFrames); 

    m_gpuUpdated = true;
}

Sample SonicSculpturePiece::sampleAudio(int windowIndex, float alpha) {
	// TODO: get window size from somewhere
	int windowSize = 512;
	// TODO: Do something better than nearest neighbors sampling...
	return m_audioSamples[windowIndex*windowSize + (int)(alpha*(windowSize-1))];
}


shared_ptr<AudioSample> SonicSculpturePiece::getBaseAudioSample() const {
    // TODO: get sample rate from somewhere
    shared_ptr<AudioSample> sound = AudioSample::createEmpty(440000);
    sound->buffer.appendPOD(m_audioSamples);
    return sound;
}

shared_ptr<AudioSample> SonicSculpturePiece::getAudioSampleFromRay(const Ray& ray) {
	Vector3 zAxis = -ray.direction();
	Vector3 xAxis;
	if (abs(zAxis.dot(Vector3::unitY())) < 0.9f) {
		xAxis = zAxis.unitCross(Vector3::unitY());
	} else {
		xAxis = zAxis.unitCross(Vector3::unitX());
	}
	Vector3 yAxis = zAxis.cross(xAxis);

	CFrame rayFrame = CFrame(Matrix3::fromColumns(xAxis, yAxis, zAxis), ray.origin());
	debugPrintf("%s\n", rayFrame.toXYZYPRDegreesString().c_str());

	// TODO: get sample rate from somewhere
	shared_ptr<AudioSample> sound = AudioSample::createEmpty(440000);

	getTransformedFramesFromOriginals();
	Array<CFrame> raySpaceFrames;
	for (auto& frame : m_transformedFrames) {
		raySpaceFrames.append(rayFrame.toObjectSpace(frame));
	}

	Array<Vector2> startEndZs;
	Array<float> minZs;
	Array<float> maxZs;
	float maximumZ = -finf();
	float minimumZ = finf();
	for (int i = 0; i < raySpaceFrames.size() - 1; ++i) {
		Vector2 startEndZ = Vector2(raySpaceFrames[i].translation.z, raySpaceFrames[i+1].translation.z);
		startEndZs.append(startEndZ);
		float minZ = min(startEndZ[0], startEndZ[1]);
		float maxZ = max(startEndZ[0], startEndZ[1]);
		minZs.append(minZ);
		maxZs.append(maxZ);
		minimumZ = min(minimumZ, minZ);
		maximumZ = max(maximumZ, maxZ);
	}
	// TODO: get these quantities from single source, right now duplicating these values;
	double delta = 0.1;
	double metersPerSample = delta / 512.0;
	int numSamples = max((int)(-minimumZ / metersPerSample), 0);


	debugPrintf("MinZ, MaxZ, numSamples: %f, %f, %d\n", minimumZ, maximumZ, numSamples);

	sound->buffer.resize(numSamples);
    memset(sound->buffer.getCArray(), 0, numSamples*sizeof(Sample));
	// Rasterize the sounds!

	for (int i = 0; i < startEndZs.size(); ++i) {
		int startIndex	= max((int)ceil(-maxZs[i] / metersPerSample), 0);
		int endIndex	= min((int)floor(-minZs[i] / metersPerSample), numSamples);
		float startZ	= startEndZs[i][0];
		float endZ		= startEndZs[i][1];
		for (int j = startIndex; j < endIndex; ++j) {
			float sampleZ = -j*metersPerSample;
			float alpha = (sampleZ - startZ) / (endZ - startZ);
			if (alpha >= 0 && alpha <= 1) {
				sound->buffer[j] = sampleAudio(i, alpha);
			}
		}
	}


	return sound;
}

shared_ptr<SonicSculpturePiece> SonicSculpturePiece::create(shared_ptr<UniversalMaterial> material) {
    shared_ptr<SonicSculpturePiece> s(new SonicSculpturePiece());
    s->m_material = material;
    return s;
}

void SonicSculpturePiece::serialize(BinaryOutput& output) const {
    m_frame.serialize(output);
    G3D::serialize(m_audioSamples, output);
    G3D::serialize(m_deltas, output);
    //G3D::serialize(m_originalFrames, output);
    G3D::serialize(m_radii, output);
}

shared_ptr<SonicSculpturePiece> SonicSculpturePiece::fromBinaryInput(shared_ptr<UniversalMaterial> material, BinaryInput& input) {
    shared_ptr<SonicSculpturePiece> s;
    s->m_frame.deserialize(input);
    deserialize<float>(s->m_audioSamples, input);
    deserialize<float>(s->m_deltas, input);
    //deserialize<CFrame>(s->m_originalFrames, input);
    deserialize<float>(s->m_radii, input);
    s->m_material = material;
    s->uploadToGPU();
    return s; 
}

void SonicSculpturePiece::insert(const CFrame & frame, const float radius, const float delta, const Array<float>& newSamples ) {
    m_gpuUpdated = false;
    m_originalFrames.append(frame);
    m_radii.append(radius);
    m_deltas.append(delta);
	m_audioSamples.append(newSamples);
    uploadToGPU();
}

void SonicSculpturePiece::render(RenderDevice * rd, const LightingEnvironment& env) {
    if (m_radii.size() < 2) {
        return;
    }
    Args args;
    env.setShaderArgs(args);
    rd->setObjectToWorldMatrix(m_frame);
    args.setMacro("OPAQUE_PASS", 1);
    setShaderArgs(args);
    rd->setCullFace(CullFace::NONE);
    LAUNCH_SHADER("UniversalSurface/UniversalSurface_render.*", args);
}



void SonicSculpturePiece::setShaderArgs(Args & args) {
    if (!m_gpuUpdated) {
        uploadToGPU();
    }
    args.setUniform("boneMatrixTexture", m_boneTexture, Sampler::buffer());
    args.setAttributeArray("g3d_Vertex", m_vertices);
    args.setAttributeArray("g3d_Normal", m_normals);
    args.setAttributeArray("g3d_TexCoord0", m_texCoord);
    args.setAttributeArray("g3d_BoneWeights", m_boneWeights);
    args.setAttributeArray("g3d_BoneIndices", m_boneIndices);
    args.setIndexStream(m_indices);
    args.setMacro("NUM_BONES", m_transformedFrames.size());
    args.setMacro("NUM_LIGHTMAP_DIRECTIONS", 0);
    args.setMacro("HAS_VERTEX_COLOR", 0);
    args.setMacro("HAS_ALPHA", m_material->hasAlpha());
    args.setMacro("ALPHA_HINT", m_material->alphaHint());
    args.setMacro("HAS_TRANSMISSIVE", m_material->hasTransmissive());
    args.setMacro("HAS_EMISSIVE", false);

    m_material->setShaderArgs(args, "material.");
}
