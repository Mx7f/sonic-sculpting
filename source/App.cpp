/** \file App.cpp */
#include "App.h"

// Tells C++ to invoke command-line main() function even on OS X and Win32.
G3D_START_AT_MAIN();

int main(int argc, const char* argv[]) {
    (void)argc; (void)argv;
    GApp::Settings settings(argc, argv);
    
    // Change the window and other startup parameters by modifying the
    // settings class.  For example:
    settings.window.width       = 1280; 
    settings.window.height      = 720;
    settings.window.asynchronous = false;
    settings.window.resizable = true;
    settings.window.caption = "Sonic Sculpting";
    settings.dataDir = FileSystem::currentDirectory();
    settings.screenshotDirectory = "../journal/";

    settings.renderer.deferredShading = false;
    settings.renderer.orderIndependentTransparency = true;


    return App(settings).run();
}


App::App(const GApp::Settings& settings) : GApp(settings) {
    renderDevice->setColorClearValue(Color3::white());
}

void App::onCleanup() {
  m_rtAudio.stopStream();
  if( m_rtAudio.isStreamOpen() )
    m_rtAudio.closeStream();
}

int audioCallback( void * outputBuffer, void * inputBuffer, unsigned int numFrames,
            double streamTime, RtAudioStreamStatus status, void * data ) {
  float * output = (float *)outputBuffer;
  float * input  = (float *)inputBuffer;
  size_t numBytes = numFrames * sizeof(float);
  memset(output, 0, numBytes);
  // Copy input buffer, handle race condition by not caring about it.
  memcpy(g_currentAudioBuffer.getCArray(), input, numBytes);
  return 0;
}



void App::initializeAudio() {

  unsigned int bufferByteCount = 0;
  unsigned int bufferFrameCount = 512;
    
  // Check for audio devices
  if( m_rtAudio.getDeviceCount() < 1 ) {
    // None :(
    debugPrintf("No audio devices found!\n");
    exit( 1 );
  }
  RtAudio::DeviceInfo info;

  unsigned int devices = m_rtAudio.getDeviceCount();


  // Let RtAudio print messages to stderr.
  m_rtAudio.showWarnings( true );

  // Set input and output parameters
  RtAudio::StreamParameters iParams, oParams;
  iParams.deviceId = m_rtAudio.getDefaultInputDevice();
  iParams.nChannels = m_audioSettings.numChannels;
  iParams.firstChannel = 0;
  oParams.deviceId = m_rtAudio.getDefaultOutputDevice();
  oParams.nChannels = m_audioSettings.numChannels;
  oParams.firstChannel = 0;
    
  // Create stream options
  RtAudio::StreamOptions options;

  g_currentAudioBuffer.resize(bufferFrameCount);
  try {
    // Open a stream
    m_rtAudio.openStream( &oParams, &iParams, m_audioSettings.rtAudioFormat, m_audioSettings.sampleRate, &bufferFrameCount, &audioCallback, (void *)&bufferByteCount, &options );
  } catch( RtAudioError& e ) {
    // Failed to open stream
    std::cout << e.getMessage() << std::endl;
    exit( 1 );
  }
  g_currentAudioBuffer.resize(bufferFrameCount);
  m_rtAudio.startStream();

}

void App::onInit() {
    GApp::onInit();
    m_appMode = AppMode::DEFAULT;
    // TODO: Print instructions
    m_maxSavedTimeSlices = 512;
    initializeAudio();

    m_rawAudioTexture = Texture::createEmpty("Raw Audio Texture", g_currentAudioBuffer.size(), 1, ImageFormat::R32F());
    m_frequencyAudioTexture = Texture::createEmpty("Frequency Audio Texture", g_currentAudioBuffer.size()/2, 1, ImageFormat::RG32F());

    m_smoothedRootMeanSquare = 0.0f;

    makeGUI();
    loadScene("Sculpting");
}

void App::makeGUI() {
    // Initialize the developer HUD (using the existing scene)
    createDeveloperHUD();
    debugWindow->setVisible(false);
    developerWindow->videoRecordDialog->setEnabled(true);
    debugPane->pack();


    debugWindow->pack();
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
    developerWindow->cameraControlWindow->setVisible(false);
    developerWindow->sceneEditorWindow->setVisible(false);
    developerWindow->setVisible(false);
    showRenderingStats = false;
    developerWindow->cameraControlWindow->moveTo(Point2(developerWindow->cameraControlWindow->rect().x0(), 0));


}


bool App::onEvent(const GEvent& e) {
    if (GApp::onEvent(e)) {
        return true;
    }
    // If you need to track individual UI events, manage them here.
    // Return true if you want to prevent other parts of the system
    // from observing this specific event.
    //
    // For example,
    // if ((e.type == GEventType::GUI_ACTION) && (e.gui.control == m_button)) { ... return true;}
    // if ((e.type == GEventType::KEY_DOWN) && (e.key.keysym.sym == GKey::TAB)) { ... return true; }

    return false;
}

void App::setAudioShaderArgs(Args& args) {
    m_rawAudioTexture->setShaderArgs(args, "rawAudio_", Sampler::video());
    m_frequencyAudioTexture->setShaderArgs(args, "frequencyAudio_", Sampler::video());
    m_fastMovingAverage.gpuData->setShaderArgs(args, "fastEWMAfreq_", Sampler::video());
    m_slowMovingAverage.gpuData->setShaderArgs(args, "slowEWMAfreq_", Sampler::video());
    m_glacialMovingAverage.gpuData->setShaderArgs(args, "glacialEWMAfreq_", Sampler::video());
}

void App::updateAudioData() {
    int sampleCount = g_currentAudioBuffer.size();
    int freqCount = sampleCount / 2;
    m_cpuRawAudioData.appendPOD(g_currentAudioBuffer);

    float sumSquare = 0.0f;
    for (int i = m_cpuRawAudioData.size() - sampleCount; i < m_cpuRawAudioData.size(); ++i) {
        sumSquare += square(m_cpuRawAudioData[i]);
    }
    float rms = sqrt(sumSquare / sampleCount);
    m_smoothedRootMeanSquare = lerp(rms, m_smoothedRootMeanSquare, 0.8f);

    int numStoredTimeSlices = m_cpuRawAudioData.size() / sampleCount;
    shared_ptr<CPUPixelTransferBuffer> ptb = CPUPixelTransferBuffer::fromData(sampleCount, numStoredTimeSlices, ImageFormat::R32F(), m_cpuRawAudioData.getCArray());
    m_rawAudioTexture->resize(sampleCount, numStoredTimeSlices);
    m_rawAudioTexture->update(ptb);

    Array<complex> frequency;
    frequency.resize(freqCount);
    float* currentRawAudioDataPtr = m_cpuRawAudioData.getCArray() + (m_cpuRawAudioData.size() - sampleCount);
    memcpy(frequency.getCArray(), currentRawAudioDataPtr, sizeof(float)*sampleCount);
    rfft((float*)frequency.getCArray(), frequency.size(), FFT_FORWARD);
    m_cpuFrequencyAudioData.appendPOD(frequency);

    Array<float> frequencyMagnitude;
    for (complex c : frequency) {
        frequencyMagnitude.append(cmp_abs(c));
    }
    if (isNull(m_fastMovingAverage.gpuData)) {
        m_fastMovingAverage.init(0.6, frequencyMagnitude, "Fast Freq EWMA");
        m_slowMovingAverage.init(0.85, frequencyMagnitude, "Slow Freq EWMA");
        m_glacialMovingAverage.init(0.95, frequencyMagnitude, "Glacial Freq EWMA");
    } else {
        m_fastMovingAverage.update(frequencyMagnitude);
        m_slowMovingAverage.update(frequencyMagnitude);
        m_glacialMovingAverage.update(frequencyMagnitude);
    }

    shared_ptr<CPUPixelTransferBuffer> freqPTB = CPUPixelTransferBuffer::fromData(freqCount, numStoredTimeSlices, ImageFormat::RG32F(), m_cpuFrequencyAudioData.getCArray());

    m_frequencyAudioTexture->resize(freqCount, numStoredTimeSlices);
    m_frequencyAudioTexture->update(freqPTB);

    if (numStoredTimeSlices == m_maxSavedTimeSlices) {
        int newTotalSampleCount = sampleCount*(numStoredTimeSlices - 1);
        int newTotalFrequencyCount = (sampleCount / 2)*(numStoredTimeSlices - 1);
        for (int i = 0; i < newTotalSampleCount; ++i) {
            m_cpuRawAudioData[i] = m_cpuRawAudioData[i + sampleCount];
        }
        m_cpuRawAudioData.resize(newTotalSampleCount);
        for (int i = 0; i < newTotalFrequencyCount; ++i) {
            m_cpuFrequencyAudioData[i] = m_cpuFrequencyAudioData[i + freqCount];
        }
        m_cpuFrequencyAudioData.resize(newTotalFrequencyCount);
    }
}

void App::onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& allSurfaces) {
    // This implementation is equivalent to the default GApp's. It is repeated here to make it
    // easy to modify rendering. If you don't require custom rendering, just delete this
    // method from your application and rely on the base class.

    if (!scene()) {
        return;
    }

    m_gbuffer->setSpecification(m_gbufferSpecification);
    m_gbuffer->resize(m_framebuffer->width(), m_framebuffer->height());
    m_gbuffer->prepare(rd, activeCamera(), 0, -(float)previousSimTimeStep(), m_settings.depthGuardBandThickness, m_settings.colorGuardBandThickness);

    m_renderer->render(rd, m_framebuffer, m_depthPeelFramebuffer, scene()->lightingEnvironment(), m_gbuffer, allSurfaces);

    // Debug visualizations and post-process effects
    rd->pushState(m_framebuffer); {
        // Call to make the App show the output of debugDraw(...)
        rd->setProjectionAndCameraMatrix(activeCamera()->projection(), activeCamera()->frame());

        for (auto& piece : m_sonicSculpturePieces) {
            piece->render(rd, scene()->lightingEnvironment());
        }
        if (notNull(m_currentSonicSculpturePiece)) {
            m_currentSonicSculpturePiece->render(rd, scene()->lightingEnvironment());
        }
        


        drawDebugShapes();
        const shared_ptr<Entity>& selectedEntity = (notNull(developerWindow) && notNull(developerWindow->sceneEditorWindow)) ? developerWindow->sceneEditorWindow->selectedEntity() : shared_ptr<Entity>();
        scene()->visualize(rd, selectedEntity, allSurfaces, sceneVisualizationSettings());

        // Post-process special effects
        m_depthOfField->apply(rd, m_framebuffer->texture(0), m_framebuffer->texture(Framebuffer::DEPTH), activeCamera(), m_settings.depthGuardBandThickness - m_settings.colorGuardBandThickness);

        m_motionBlur->apply(rd, m_framebuffer->texture(0), m_gbuffer->texture(GBuffer::Field::SS_EXPRESSIVE_MOTION),
            m_framebuffer->texture(Framebuffer::DEPTH), activeCamera(),
            m_settings.depthGuardBandThickness - m_settings.colorGuardBandThickness);
    } rd->popState();

    if ((submitToDisplayMode() == SubmitToDisplayMode::MAXIMIZE_THROUGHPUT) && (!renderDevice->swapBuffersAutomatically())) {
        // We're about to render to the actual back buffer, so swap the buffers now.
        // This call also allows the screenshot and video recording to capture the
        // previous frame just before it is displayed.
        swapBuffers();
    }

    // Clear the entire screen (needed even though we'll render over it, since
    // AFR uses clear() to detect that the buffer is not re-used.)
    rd->clear();

    // Perform gamma correction, bloom, and SSAA, and write to the native window frame buffer
    m_film->exposeAndRender(rd, activeCamera()->filmSettings(), m_framebuffer->texture(0));
}

void App::onUserInput(UserInput* ui) {
    GApp::onUserInput(ui);
    if (ui->keyDown(GKey::SPACE)) {
        m_appMode = AppMode::MAKING_SCULPTURE;
    } else {
        m_appMode = AppMode::DEFAULT;
    }

    // Add key handling here based on the keys currently held or
    // ones that changed in the last frame.
}

void App::onGraphics2D(RenderDevice* rd, Array<Surface2D::Ref>& posed2D) {
    // Render 2D objects like Widgets.  These do not receive tone mapping or gamma correction
    Surface2D::sortAndRender(rd, posed2D);
}

void App::onSimulation(RealTime rdt, SimTime sdt, SimTime idt) {
    GApp::onSimulation(rdt, sdt, idt);
    updateAudioData();  
    float delta = 0.1f;
    float radius = m_smoothedRootMeanSquare * 1.0f;
    CFrame frame = activeCamera()->frame();
    frame.translation += activeCamera()->frame().lookVector() * 1.1f;
    if (m_appMode == AppMode::DEFAULT) {
        if (notNull(m_currentSonicSculpturePiece)) {
            if (m_currentSonicSculpturePiece->size() > 0) {
                m_currentSonicSculpturePiece->insert(frame, 0.0f, delta);
            }
            m_sonicSculpturePieces.append(m_currentSonicSculpturePiece);
            m_currentSonicSculpturePiece = shared_ptr<SonicSculpturePiece>();
        }
    } else if (m_appMode == AppMode::MAKING_SCULPTURE) {
        if (isNull(m_currentSonicSculpturePiece)) {
            m_currentSonicSculpturePiece = SonicSculpturePiece::create(UniversalMaterial::createDiffuse(Color3(0.9f)));
        }
        m_currentSonicSculpturePiece->insert(frame, radius, delta);
    }

    // Example GUI dynamic layout code.  Resize the debugWindow to fill
    // the screen horizontally.
    debugWindow->setRect(Rect2D::xywh(0, 0, (float)window()->width(), debugWindow->rect().height()));
}
