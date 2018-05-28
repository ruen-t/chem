class App {
  constructor($state) {
    console.log("App controller")
    var self = this;
    self.message = "first page"
    $state.go('reaction');
  }
}

angular
  .module('app')
  .component('app', {
    templateUrl: 'app/containers/App.html',
    controller: App,
    controllerAs: "vm"
  })
 
