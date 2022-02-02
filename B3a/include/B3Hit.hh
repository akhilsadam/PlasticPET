//============ header file =====================

#ifndef B3Hit_h
#define B3Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class B3Hit : public G4VHit
{
  public:
      B3Hit();
      ~B3Hit();
      B3Hit(const B3Hit &right);
      const B3Hit& operator=(const B3Hit &right);
      int operator==(const B3Hit &right) const;

      inline void * operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw() const;
      void Print() const;

  private:
      G4double Edep;
      G4ThreeVector pos;

  public:
      inline void SetEdep(G4double edep)
      { Edep = edep; }
      inline G4double GetEdep() const
      { return Edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos() const
      { return pos; }
};

typedef G4THitsCollection<B3Hit> B3HitsCollection;

extern G4ThreadLocal G4Allocator<B3Hit>* B3HitAllocator;

inline void* B3Hit::operator new(size_t)
{
  if(!B3HitAllocator) B3HitAllocator = new G4Allocator<B3Hit>;
  return (void *) B3HitAllocator->MallocSingle();
}

inline void B3Hit::operator delete(void *aHit)
{
  B3HitAllocator->FreeSingle((B3Hit*) aHit);
}

#endif
