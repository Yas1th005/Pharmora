import '@google/model-viewer';
import model from "../imgaes/glucosemolecule.glb"

export default function LandModel() {
  return (
    <div style={{ width: '500px', height: '600px', display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
      <model-viewer
        src={model}
        alt="A 3D model of something"
        auto-rotate
        camera-controls
        style={{ width: '100%', height: '100%' }}
      ></model-viewer>
    </div>
  );
}
