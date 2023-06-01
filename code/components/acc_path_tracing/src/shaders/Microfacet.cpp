#include "shaders/Microfacet.hpp"
#include "samplers/SamplerInstance.hpp"
#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/epsilon.hpp"
#include "glm/gtc/constants.hpp"
#include "glm/gtx/compatibility.hpp"
#include "Onb.hpp"
#include "algorithm"
const float F0=0.04;
const float metalness = 0.2;
const float roughness = 0.2;
namespace AccPathTracer
{
    struct brdfData {
        Vec3 specularF0;
        Vec3 diffuseReflectance;
        Vec3 base_color;
        float alpha;        //< linear roughness - often 'alpha' in specular BRDF equations
        float alphaSquared; //< alpha squared - pre-calculated value commonly used in BRDF equations
        // Roughnesses
        float roughness;    //< perceptively linear roughness (artist's input)
        float metalness;

        // Commonly used terms for BRDF evaluation
        Vec3 F; //< Fresnel term

        // Vectors
        Vec3 V; //< Direction to viewer (or opposite direction of incident ray)
        Vec3 N; //< Shading normal
        Vec3 H; //< Half vector (microfacet normal)
        Vec3 L; //< Direction to light (or direction of reflecting ray)

        float NdotL;
        float NdotV;

        float LdotH;
        float NdotH;
        float VdotH;

        // True when V/L is backfacing wrt. shading normal N
        bool Vbackfacing;
        bool Lbackfacing;
    };
    inline float saturate(float x) { return clamp(x, 0.0f, 1.0f); };

    inline Vec3 evalFresnel(Vec3 f0, float LdotH) {
        return f0 + (Vec3{1.0f,1.0f,1.0f} - f0) * pow(1.0f - max(0.00001f,LdotH), 5.0f);
    };

    inline float GGX_D(float alphaSquared, float NdotH) {
        float b = ((alphaSquared - 1.0f) * NdotH * NdotH + 1.0f);
        return alphaSquared / (PI * b * b);
    };

    inline float Smith_G_a(float NdotS,float alpha) {
        return NdotS / (max(0.00001f, alpha) * sqrt(1.0f - min(0.99999f, NdotS * NdotS)));
    };

    // Lambda function for Smith G term derived for GGX distribution
    inline float Smith_G_Lambda(float a) {
        return (-1.0f + sqrt(1.0f + (1.0f / (a * a)))) * 0.5f;
    };

    inline float Smith_G2(float alpha, float NdotL, float NdotV) {
        float aL = Smith_G_a(NdotL,alpha);
        //cout << "al:" << aL << endl;
        float aV = Smith_G_a(NdotV,alpha);
        //cout << "aV:" << aV << endl;
        //return 1.0f / (1.0f + Smith_G_Lambda(aL) + Smith_G_Lambda(aV));
        return aL * aV;
    };

    inline Vec3 evalMicrofacet(brdfData data) {
        float D = GGX_D(data.alphaSquared, data.NdotH);
        float G2 = Smith_G2(data.alpha, data.NdotL, data.NdotV);

        return data.F * (G2 * D * data.NdotL);
    };

    brdfData buildData(const Ray& ray, const Vec3& hitPoint, const Vec3& normal,RenderOption renderoption,Vec3 albedo, Vec3 direction) {
        brdfData data;
        //data.roughness = renderoption.roughness;
        //data.metalness = renderoption.metalness;
        data.roughness = roughness;
        data.metalness = metalness;
        //data.alpha = data.roughness * data.roughness;
        data.alpha = 0.5f;
        data.alphaSquared = data.alpha * data.alpha;
        //data.base_color = albedo;
        data.base_color= Vec3(0.9f, 0.4f, 0.4f);
        data.diffuseReflectance = data.base_color * (1.0f - data.metalness);
        //通过插值获取F0
        //glm::vec3 result = glm::lerp(glm::lowp_vec3(a), glm::lowp_vec3(b), glm::lowp_vec3(c));

        //data.specularF0 = (1 - data.metalness) * Vec3{ F0, F0, F0 } + data.metalness * data.base_color;
        data.specularF0 = Vec3(0.31f, 0.31f, 0.31f);

        //data.specularF0 = glm::lerp(Vec3{ F0, F0, F0 }, data.base_color, data.metalness);
        data.V = Vec3{ -(ray.direction.x),-(ray.direction.y),-(ray.direction.z) };
        data.L = direction;
        //查看随机生成的方向
        //cout << "L(" << L.x << "," << L.y << "," << L.z << ")" << endl;
        data.N = normal;
        data.H = glm::normalize(data.V+data.L);
        data.NdotL = glm::dot(data.N, data.L);
        data.NdotV = glm::dot(data.N, data.V);
        //cout << "NdotL:" << NdotL > 0 ? "true" : "false" << endl;
        //cout << "NdotV:" << NdotV > 0 ? "true" : "false" << endl;
        //Vbackfacing = (NdotV <= 0.0f);
        //Lbackfacing = (NdotL <= 0.0f);
        
        //需要确保点乘在0~1之间
        data.NdotL = glm::min(glm::max(0.0f, data.NdotL), 1.0f);
        data.NdotV = glm::min(glm::max(0.0f, data.NdotV), 1.0f);

        data.LdotH = saturate(glm::dot(data.L, data.H));
        data.NdotH = saturate(glm::dot(data.N, data.H));
        data.VdotH = saturate(glm::dot(data.V, data.H));

        //cout << "LdotH:" << glm::dot(data.V, data.H) << endl;
        data.F = evalFresnel(data.specularF0, data.VdotH);

        return data;
    }
    Microfacet::Microfacet(Material& material, vector<Texture>& textures, RenderOption& renderopt)
        : Shader(material, textures)
    {
        auto diffuseColor = material.getProperty<Property::Wrapper::RGBType>("diffuseColor");
        if (diffuseColor) albedo = (*diffuseColor).value;
        else albedo = { 1, 1, 1 };
        renderoption = renderopt;
    }
    
    Scattered Microfacet::shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal)const {
        Vec3 origin = hitPoint;
        Vec3 random = defaultSamplerInstance<HemiSphere>().sample3d();

        Onb onb{ normal };
        Vec3 direction = glm::normalize(onb.local(random));
        //cout<<"ditection:(" << direction.x << "," << direction.y << "," << direction.z << ")" << endl;
        brdfData data=buildData(ray,hitPoint,normal,renderoption,albedo,direction);

        //if (Vbackfacing || Lbackfacing) return { Ray{origin,direction},Vec3{0},Vec3{0} }
        //cout<<"F:"<< data.F.x << "," << data.F.y << "," << data.F.z << ")" << endl;

        Vec3 specular = evalMicrofacet(data);
        //Vec3 diffuse = (Vec3{ 1.0f, 1.0f, 1.0f } - data.F) * data.diffuseReflectance * data.NdotL / PI;
        Vec3 diffuse = data.base_color;
        float pdf = 1 / (2 * PI);
        return {
            Ray{origin, direction},
            diffuse,
            specular,
            pdf
        };
    }

}