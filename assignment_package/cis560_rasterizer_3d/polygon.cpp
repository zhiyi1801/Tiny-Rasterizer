#include "polygon.h"
#include <glm/gtx/transform.hpp>

void Polygon::Triangulate()
{
    //TODO: Populate list of triangles
    for(unsigned int i = 0; i < this->m_verts.size() - 2; ++i)
    {
        this->m_tris.push_back(Triangle{0,i+1,i+2});
    }
}

glm::vec3 GetImageColor(const glm::vec2 &uv_coord, const QImage* const image)
{
    if(image)
    {
        int X = glm::min(image->width() * uv_coord.x, image->width() - 1.0f);
        int Y = glm::min(image->height() * (1.0f - uv_coord.y), image->height() - 1.0f);
        QColor color = image->pixel(X, Y);
        return glm::vec3(color.red(), color.green(), color.blue());
    }
    return glm::vec3(255.f, 255.f, 255.f);
}

glm::vec3 GetBarycentric(const Vertex &v1, const Vertex &v2, const Vertex &v3, const glm::vec3 &p)
{
    glm::vec3 lAB(v1.m_pos - v2.m_pos);
    glm::vec3 lAC(v1.m_pos - v3.m_pos);
    glm::vec3 lPA(p - glm::vec3(v1.m_pos));
    glm::vec3 vecx(lAB.x, lAC.x, lPA.x);
    glm::vec3 vecy(lAB.y, lAC.y, lPA.y);
    glm::vec3 uv1 = glm::cross(vecx, vecy);
    //vertex are all int
    if(std::abs(uv1.z) < 1) return glm::vec3(-1,-1,1);
    else
    {
        uv1 /= uv1.z;
        return glm::vec3(1-uv1.x-uv1.y, uv1.x, uv1.y);
    }
}

glm::vec3 GetBarycentric(const glm::vec4 &v1, const glm::vec4 &v2, const glm::vec4 &v3, const glm::vec3 &p)
{
    glm::vec3 lAB(v1 - v2);
    glm::vec3 lAC(v1 - v3);
    glm::vec3 lPA(p - glm::vec3(v1));
    glm::vec3 vecx(lAB.x, lAC.x, lPA.x);
    glm::vec3 vecy(lAB.y, lAC.y, lPA.y);
    glm::vec3 uv1 = glm::cross(vecx, vecy);
    //vertex are all int
    if(std::abs(uv1.z) < 1) return glm::vec3(-1,-1,1);
    else
    {
        uv1 /= uv1.z;
        return glm::vec3(1-uv1.x-uv1.y, uv1.x, uv1.y);
    }
}

glm::vec3 GetBarycentric(const glm::vec4 &v1, const glm::vec4 &v2, const glm::vec4 &v3, const glm::vec4 &p)
{
    glm::vec3 lAB(v1 - v2);
    glm::vec3 lAC(v1 - v3);
    glm::vec3 lPA(p - v1);
    glm::vec3 vecx(lAB.x, lAC.x, lPA.x);
    glm::vec3 vecy(lAB.y, lAC.y, lPA.y);
    glm::vec3 uv1 = glm::cross(vecx, vecy);
    //vertex are all int
    if(std::abs(uv1.z) < 1) return glm::vec3(-1,-1,1);
    else
    {
        uv1 /= uv1.z;
        return glm::vec3(1-uv1.x-uv1.y, uv1.x, uv1.y);
    }
}

// Creates a polygon from the input list of vertex positions and colors
Polygon::Polygon(const QString& name, const std::vector<glm::vec4>& pos, const std::vector<glm::vec3>& col)
    : m_tris(), m_verts(), m_name(name), mp_texture(nullptr), mp_normalMap(nullptr)
{
    for(unsigned int i = 0; i < pos.size(); i++)
    {
        m_verts.push_back(Vertex(pos[i], col[i], glm::vec4(), glm::vec2()));
    }
    Triangulate();
}

// Creates a regular polygon with a number of sides indicated by the "sides" input integer.
// All of its vertices are of color "color", and the polygon is centered at "pos".
// It is rotated about its center by "rot" degrees, and is scaled from its center by "scale" units
Polygon::Polygon(const QString& name, int sides, glm::vec3 color, glm::vec4 pos, float rot, glm::vec4 scale)
    : m_tris(), m_verts(), m_name(name), mp_texture(nullptr), mp_normalMap(nullptr)
{
    glm::vec4 v(0.f, 1.f, 0.f, 1.f);
    float angle = 360.f / sides;
    for(int i = 0; i < sides; i++)
    {
        glm::vec4 vert_pos = glm::translate(glm::vec3(pos.x, pos.y, pos.z))
                           * glm::rotate(rot, glm::vec3(0.f, 0.f, 1.f))
                           * glm::scale(glm::vec3(scale.x, scale.y, scale.z))
                           * glm::rotate(i * angle, glm::vec3(0.f, 0.f, 1.f))
                           * v;
        m_verts.push_back(Vertex(vert_pos, color, glm::vec4(), glm::vec2()));
    }

    Triangulate();
}

Polygon::Polygon(const QString &name)
    : m_tris(), m_verts(), m_name(name), mp_texture(nullptr), mp_normalMap(nullptr)
{}

Polygon::Polygon()
    : m_tris(), m_verts(), m_name("Polygon"), mp_texture(nullptr), mp_normalMap(nullptr)
{}

Polygon::Polygon(const Polygon& p)
    : m_tris(p.m_tris), m_verts(p.m_verts), m_name(p.m_name), mp_texture(nullptr), mp_normalMap(nullptr)
{
    if(p.mp_texture != nullptr)
    {
        mp_texture = new QImage(*p.mp_texture);
    }
    if(p.mp_normalMap != nullptr)
    {
        mp_normalMap = new QImage(*p.mp_normalMap);
    }
}

Polygon::~Polygon()
{
    delete mp_texture;
}

void Polygon::SetTexture(QImage* i)
{
    mp_texture = i;
}

void Polygon::SetNormalMap(QImage* i)
{
    mp_normalMap = i;
}

void Polygon::AddTriangle(const Triangle& t)
{
    m_tris.push_back(t);
}

void Polygon::AddVertex(const Vertex& v)
{
    m_verts.push_back(v);
}

void Polygon::ClearTriangles()
{
    m_tris.clear();
}

Triangle& Polygon::TriAt(unsigned int i)
{
    return m_tris[i];
}

Triangle Polygon::TriAt(unsigned int i) const
{
    return m_tris[i];
}

Vertex &Polygon::VertAt(unsigned int i)
{
    return m_verts[i];
}

Vertex Polygon::VertAt(unsigned int i) const
{
    return m_verts[i];
}

LineSegment2D::LineSegment2D(const glm::vec3& p1, const glm::vec3& p2):pStart(p1.y<=p2.y?p1:p2),pEnd(p1.y<=p2.y?p2:p1),dy(p2.y-p1.y),dx(p2.x-p1.x){}

LineSegment2D::LineSegment2D(const glm::vec4& p1, const glm::vec4& p2):pStart(p1.y<=p2.y?p1:p2),pEnd(p1.y<=p2.y?p2:p1),dy(p2.y-p1.y),dx(p2.x-p1.x){}

LineSegment2D::LineSegment2D():pStart(glm::vec3(0)),pEnd(glm::vec3(0)),dy(0),dx(0){}

LineSegment2D::LineSegment2D(const LineSegment2D &l):pStart(l.pStart), pEnd(l.pEnd), dy(l.dy), dx(l.dx){};


bool LineSegment2D::getIntersection(int y, float* x)
{
    if(y<pStart.y || y > pEnd.y) return false;
    else if(dx == 0)
    {
        *x = pStart.x;
        return true;
    }
    else
    {
        *x = pStart.x + (y - pStart.y) * dx/dy;
        return true;
    }
}
