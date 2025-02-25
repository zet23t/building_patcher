#include <raymath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "game.h"
#include "RuntimeContext.h"
#include "shared/ui/dusk-gui.h"
#include "PiecePatcher.h"

#ifndef COMPONENT_IMPLEMENTATION
#include "shared/scene_graph/components/CameraComponent.h"
#endif

#define COMPONENT(t) void t##_register();
#include "component_list.h"
#undef COMPONENT

ComponentIdMap _componentIdMap;
RuntimeContext *_runtimeContext;
ResourceManager *_resourceManager;

static SceneGraph *_scene;

int game_init(RuntimeContext *runtimeContext)
{
    _scene = SceneGraph_create();
    _runtimeContext = runtimeContext;
    _resourceManager = &runtimeContext->resourceManager;

    PiecePatcherInit();
    DuskGui_init();
    Font font = ResourceManager_loadFont(&runtimeContext->resourceManager, "assets/myfont-regular.png");
    DuskGui_setDefaultFont(font, font.baseSize, 1);

    // Register components to the scene graph
#define COMPONENT(t) t##_register(_scene);
#include "component_list.h"
#undef COMPONENT

    SceneObjectId cameraId = SceneGraph_createObject(_scene, "Camera");
    SceneComponentId cameraComponentId = SceneGraph_addComponent(_scene, cameraId, _componentIdMap.CameraComponentId,
        &(CameraComponent) {
            .camera.near = 0.1f,
            .camera.far = 100.0f,
            .camera.fovy = 15.0f,
            .camera.projection = CAMERA_ORTHOGRAPHIC,
            .targetDistance = 10.0f,
            .targetPoint = (Vector3){0.0f, 0.0f, 0.0f},
            .trackTarget = 1,
        });
    SceneGraph_setLocalRotation(_scene, cameraId, (Vector3){15.0f, 170.0f, 0.0f});

    // SceneObjectId playerObjectId = SceneGraph_createObject(_scene, "Cube");
    // SceneGraph_addComponent(_scene, playerObjectId, _componentIdMap.PrimitiveRendererComponentId,
    //     &(PrimitiveRendererComponent) {
    //         .primitiveType = PRIMITIVE_TYPE_CUBE,
    //         .color = (Color){255, 0, 0, 255},
    //         .size = (Vector3){1.0f, 1.0f, 1.0f},
    //     });

    SceneObjectId systemsId = SceneGraph_createObject(_scene, "Systems");
    SceneGraph_addComponent(_scene, systemsId, _componentIdMap.ChunkSystemComponentId,
        &(ChunkSystemComponent) {
            .chunkLength = 4,
            .cameraComponentId = cameraComponentId,
        });
    // SceneGraph_addComponent(_scene, systemsId, _componentIdMap.GameplaySystemComponentId,
    //     &(GameplaySystemComponent) {
    //         .playerObjectId = playerObjectId,
    //         .cameraComponentId = cameraComponentId,
    //     });

    return 1;
}

void game_draw()
{
    Camera3D camera = CameraComponent_getCamera3D(_scene, (SceneObjectId){0});
    BeginMode3D(camera);
    // SceneGraph_draw(_scene, camera, NULL);
    // SceneGraph_sequentialDraw(_scene, camera, NULL);
    EndMode3D();

    PiecePatcherDrawDebug();
    

    // int screenWidth = GetScreenWidth();
    // int screenHeight = GetScreenHeight();
    // int panelHeight = 150;
    // int panelWidth = 350;
    // DuskGuiParamsEntryId panel = DuskGui_beginPanel(
    //     (DuskGuiParams){
    //         .rayCastTarget = 1,
    //         .bounds = (Rectangle){
    //             .x = (screenWidth - panelWidth) * 0.5f, 
    //             .y = (screenHeight - panelHeight) - 10, 
    //             .width = panelWidth, .height = panelHeight}});

    // static DuskGuiStyleGroup group = {0};
    // group = *DuskGui_getStyleGroup(DUSKGUI_STYLE_LABEL);
    // group.fallbackStyle.textAlignment.y = 0.0f;

    // DuskGui_label((DuskGuiParams){
    //     .text = "Hello, World!\n"
    //         "This UI is made with the [color=547f]Dusk GUI system[/color].\n"
    //         "You can recompile and reload with [color=00ff]F8[/color].\n"
    //         "The [color=830f]editor[/color] UI is currently only a debugging tool.\n",
    //     .bounds = (Rectangle){10.0f, 10.0f, panelWidth, panelHeight},
    //     .styleGroup = &group
    // });

    // DuskGui_endPanel(panel);
    
    Editor_draw(&_runtimeContext->editorState, _scene);
    Editor_drawControls(&_runtimeContext->editorState, _scene);

    DuskGui_finalize();
}

void game_update(float dt)
{
    if (_runtimeContext->editorState.paused) return;
    _runtimeContext->editorState.gameTime += dt;
    _runtimeContext->editorState.frameCount += 1;

    SceneGraph_updateTick(_scene, dt);
}