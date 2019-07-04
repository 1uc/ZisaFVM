#ifndef ZISA_FACE_HPP_JIDEW
#define ZISA_FACE_HPP_JIDEW

namespace zisa {
class Face {
public:
  DenormalizedRule qr;

  Face() = default;
  explicit Face(DenormalizedRule qr) : qr(std::move(qr)) {}
};

inline bool operator==(const Face &a, const Face &b) { return a.qr == b.qr; }
inline bool operator!=(const Face &a, const Face &b) { return !(a == b); }

inline double volume(const Face &face) { return volume(face.qr); }

}

#endif // ZISA_FACE_HPP
