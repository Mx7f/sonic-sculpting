#include "SoundInstance.h"

void SoundInstance::initialize() {
    int numSamplesPerRow = 512; // TODO: Pull from somewhere
    int numAudioRows = (audioSample->buffer.size() + (numSamplesPerRow - 1)) / numSamplesPerRow;
    Array<float> paddedSamples;
    paddedSamples.appendPOD(audioSample->buffer);
    int numSamples = audioSample->buffer.size();
    while (paddedSamples.size() < numSamplesPerRow * numAudioRows) {
        paddedSamples.append(0.0f);
    }
    shared_ptr<CPUPixelTransferBuffer> ptb = CPUPixelTransferBuffer::fromData(numSamplesPerRow, numAudioRows, ImageFormat::R32F(), paddedSamples.getCArray());
    m_rawAudioTexture = Texture::createEmpty("SoundInstance Raw Audio " + name, numSamplesPerRow, numAudioRows, ImageFormat::R32F());
    m_rawAudioTexture->update(ptb);
    int width = 1024;
    int height = 400;
    m_displayTexture = Texture::createEmpty("Sound Instance " + name, width, height);
    m_framebuffer = Framebuffer::create(m_displayTexture);
    RenderDevice* rd = RenderDevice::current;
    rd->push2D(m_framebuffer); {
        rd->setColorClearValue(Color4(0.0, 0.0, 0.0, 0.5));
        rd->clear();
        rd->setBlendFunc(RenderDevice::BLEND_ONE, RenderDevice::BLEND_ONE, RenderDevice::BLENDEQ_ADD, RenderDevice::BLENDEQ_SAME_AS_RGB, Framebuffer::COLOR0);
        Color3 color = Color3::white()*0.8f;
        float yOffset = 0.0f;
        Args args;
        args.setUniform("color", color);
        args.setUniform("waveformScale", Vector2(width* 0.5f, height* 0.45f));
        args.setUniform("offset", Vector3(width * 0.5f, height * 0.5f, 0.0f));
        args.setUniform("numSamples", numSamples);
        m_rawAudioTexture->setShaderArgs(args, "rawAudio_", Sampler::buffer());
        args.setPrimitiveType(PrimitiveType::LINE_STRIP);
        args.setNumIndices(numSamples);
        LAUNCH_SHADER("visualizeLines.*", args);
    } rd->pop2D();
}
