#include "server/Server.hpp"
#include "scene/Scene.hpp"
#include "component/RenderComponent.hpp"
#include "Camera.hpp"
#include "Metropolis.hpp"

using namespace std;
using namespace NRenderer;

namespace Metropolis
{

    class Adapter : public RenderComponent
    {
        void render(SharedScene spScene) {
            MetropolisRenderer renderer{ spScene };
            auto renderResult = renderer.render();
            auto [pixels, width, height] = renderResult;
            getServer().screen.set(pixels, width, height);
            renderer.release(renderResult);

            // logger
            //getServer().logger.log("common...");
            //getServer().logger.success("success...");
            //getServer().logger.warning("warning...");
            //getServer().logger.error("error...");
        }
    };
}

// REGISTER_RENDERER(Name, Description, Class)
REGISTER_RENDERER(MetropolisLightTransport, "MLT rendering algorithm.", Metropolis::Adapter);